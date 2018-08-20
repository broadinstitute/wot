import numexpr as ne
import numpy as np
from scipy.optimize import fmin_tnc
from scipy.stats import entropy
from sklearn.metrics.pairwise import pairwise_distances


def unfold(tensor, mode):
    xm = np.rollaxis(tensor, mode, 0).reshape(tensor.shape[mode], -1).T
    return xm


def fold(xm, mode, new_shape):
    tmp_shp = tuple([new_shape[mode]] + [new_shape[i] for i in range(len(new_shape)) if i != mode])
    xt = xm.T.reshape(tmp_shp)
    xt = np.rollaxis(xt, 0, mode + 1)
    return xt


def ttm(tensor, matrix, mode):
    xm = unfold(tensor, mode).dot(matrix)
    new_shape = list(tensor.shape)
    new_shape[mode] = matrix.shape[1]
    return fold(xm, mode, tuple(new_shape))


class SparseOptimization(object):
    def __init__(self, threads):
        self.totally_radical = True
        _ = ne.set_num_threads(threads)
        self.threads = threads
        self.Z_loss_series = []
        self.U_loss_series = []
        self.Z_grad_series = []
        self.U_grad_series = []

    def set_fa(self, k, b, y0, x0):
        def gen_log(x, k=k, b=b, y0=y0, x0=x0):
            return ne.evaluate("1./(1 + exp(-(x-x0)*b)*(k/y0-1))*k")

        def gen_log_grad(x, k=k, b=b, y0=y0, x0=x0):
            fx = gen_log(x, k, b, y0, x0)
            return ne.evaluate("(1 - fx/k)*fx*b")

        self.fa = gen_log
        self.fa_grad = gen_log_grad
        self.k, self.b, self.y0, self.x0 = k, b, y0, x0

    def nonlinear_proxGrad(self, params, maxItr=100, tol=1e-5, step_factor=2, acceptance_const=0.5, crit_size=1,
                           eps=1e-3, with_prints=False, fa_update_freq=10, nonneg=False, forceNorm=False, printFreq=10):
        for itr in range(maxItr):
            params = self.params_series[-1]
            grad = self.get_grad(params)
            self.grad_series = self.grad_series[-1:] + [grad]
            if forceNorm:
                step = self.get_step(min_step=1e-2)
            else:
                step = self.get_step(min_step=1.)
            self.loss_series.append(self.simple_loss(self.params_series[-1]))
            self.loss_series = self.loss_series[-crit_size:]
            while True:
                update = params - grad / step
                p = self.proximal_optimum(update, self.lda1 / step, nonneg=nonneg, forceNorm=forceNorm)
                # acceptance criterion
                cond = 0
                ac = acceptance_const * step / 2 * np.linalg.norm(p - params) ** 2
                for l in self.loss_series:
                    c = l - ac
                    if c > cond:
                        cond = c
                slz = self.simple_loss(p)
                # print itr,step,slz,cond,self.loss_series,ac
                if (cond == 0) or (slz <= cond * (1 + eps)):
                    self.params_series = self.params_series[-1:] + [p]
                    break
                step *= step_factor
                if step > 2 ** 100:
                    break
            if step > 2 ** 100:
                print('FAILURE: prox grad step size -> 0. Resetting to initial value...')
                params = self.params_series[0]
                if (itr < fa_update_freq) and (fa_update_freq < maxItr):
                    print('updating fa')
                    self.update_fa(params)
                break
            params = self.params_series[-1]
            if (itr > 0) and (itr % fa_update_freq == 0):
                print('updating fa')
                self.update_fa(params)
                # add method to update loss series here for new fa
                self.update_loss_series()
            if with_prints and (itr % printFreq == 0):
                self.print_performance(params, itr, self.params_series)
            if ((np.linalg.norm(self.params_series[-2] - self.params_series[-1]) / np.linalg.norm(
                    self.params_series[-2]) < tol) or (np.linalg.norm(self.params_series[-2]) == 0)) and (itr > 2):
                if (itr < fa_update_freq) and (fa_update_freq < maxItr):
                    print('updating fa')
                    self.update_fa(params)
                break
        return params

    def print_performance(self, params, itr, params_series=(1., 1.)):
        loss = 2 * self.simple_loss(params) - self.lda1 * abs(params).sum() - self.lda2 * np.linalg.norm(params) ** 2
        if self.withCoupling:
            print('itr: %d, loss: %f, params entropy: %f, params min: %f, params max: %f, params change: %f' % (
            itr, loss, np.average([np.exp(entropy(abs(z))) for z in params.T]), params.min(), params.max(),
            np.linalg.norm(params_series[-2] - params_series[-1]) / np.linalg.norm(params_series[-2])))
        else:
            fit = 1 - loss / sum([np.linalg.norm(x) ** 2 / x.shape[0] for x in self.Xg[1:] if len(x) > 0])
            print('itr: %d, fit: %f, params entropy: %f, params min: %f, params max: %f, params change: %f' % (
            itr, fit, np.average([np.exp(entropy(abs(z))) for z in params.T]), params.min(), params.max(),
            np.linalg.norm(params_series[-2] - params_series[-1]) / np.linalg.norm(params_series[-2])))

    def get_step(self, min_step=0.0001):
        if len(self.grad_series) > 1:
            d = self.params_series[-1] - self.params_series[-2]
            g = self.grad_series[-1] - self.grad_series[-2]
            return max(min_step, min((d * g).sum() / (d * d).sum(), (g * g).sum() / (d * g).sum()) / 100.)
        else:
            return min_step

    def simple_loss(self, params):
        loss = 0
        for t, lineage, xg, xr in zip(range(1, len(self.Xg)), self.Lineage, self.Xg[1:], self.Xr[:-1]):
            if len(lineage) > 0:
                Xh = self.get_Xhat(params, xr)
                if self.withCoupling:
                    xgt1 = self.Xg[t - 1]
                    x = ne.evaluate("xgt1 + Xh")
                    d = pairwise_distances(x, Y=xg, metric='sqeuclidean', n_jobs=self.threads)
                    loss += ne.evaluate("d*lineage").sum()
                else:
                    loss += np.linalg.norm(ne.evaluate("xg - Xh")) ** 2 / xg.shape[0]
        return 0.5 * loss + self.lda1 * abs(params).sum() + self.lda2 * np.linalg.norm(params) ** 2

    def get_all_Xhat(self, params):
        Xh = []
        for xr in self.Xr:
            if len(xr) > 0:
                Xh.append(self.get_Xhat(params, xr))
            else:
                Xh.append([])
        return Xh

    def proximal_optimum(self, params, delta, nonneg=False, forceNorm=False):
        if delta > 0:
            z = (params - delta * np.sign(params)) * (abs(params) > delta)
        else:
            z = params
        if nonneg:
            z[(z < 0)] = 0
        elif hasattr(self, 'prox_bounds'):
            z = np.maximum(z, self.prox_bounds[0])
            z = np.minimum(z, self.prox_bounds[1])
        if forceNorm:
            if z.shape[0] < z.shape[1]:
                z = (z.T / np.linalg.norm(z, axis=1)).T
            else:
                z = z / np.linalg.norm(z, axis=0)
            z[np.isnan(z)] = 0
        return z

    def get_W(self, params, xr, do_grad=False, func_args={}, return_linear=False, both=False):
        W = xr.dot(params)
        if return_linear:
            return W
        if both:
            Wfa = self.fa(W, **func_args)
            Wg = self.fa_grad(W, **func_args)
            return W, Wfa, Wg
        if do_grad:
            W = self.fa_grad(W, **func_args)
        else:
            W = self.fa(W, **func_args)
        return W

    def grad_Z(self, Z):
        gZ = np.zeros(Z.shape)
        for t, lineage, xg, xr in zip(range(1, len(self.Xg)), self.Lineage, self.Xg[1:], self.Xr[:-1]):
            if len(lineage) > 0:
                Wlinear, W, Wgrad = self.get_W(Z, xr, do_grad=True, both=True)
                if self.withCoupling:
                    xhu = W.dot(self.U).dot(self.U.T)
                    xgu = xg.dot(self.U.T)[:, np.newaxis] - self.Xg[t - 1].dot(self.U.T)
                    module_resid = xgu - xhu
                    module_resid = - module_resid
                    # module_resid is now (t2 cells x t1 cells x modules)
                    # now we weight by lineage info
                    module_resid = module_resid.T
                    module_resid = ne.evaluate("module_resid*lineage").T
                    module_resid = ne.evaluate("module_resid*Wgrad")
                    gZ += ttm(module_resid, xr, 1).sum(0)
                else:
                    xh = W.dot(self.U)
                    module_resid = ne.evaluate("xh - xg").dot(self.U.T)
                    module_resid = ne.evaluate("module_resid*Wgrad")
                    gZ += xr.T.dot(module_resid) / xg.shape[0]
        return gZ + self.lda2 * Z

    def grad_U(self, U):
        gU = np.zeros(U.shape)
        for t, lineage, xg, xr in zip(range(1, len(self.Xg)), self.Lineage, self.Xg[1:], self.Xr[:-1]):
            if len(lineage) > 0:
                W = self.get_W(self.Z, xr)
                xh = W.dot(U)
                if self.withCoupling:
                    for i in range(xg.shape[0]):
                        for j in range(xr.shape[0]):
                            r = xh[j] - (xg[i] - self.Xg[t - 1][j])
                            r *= lineage[j, i]
                            gU += np.multiply.outer(W[j], r)
                else:
                    gU += W.T.dot(ne.evaluate("xh - xg")) / xg.shape[0]
        return gU + self.lda2 * U

    def update_Z(self, maxItr=100, with_prints=False, fa_update_freq=10, nonneg=False, forceNorm=False):
        if with_prints:
            print('updating Z')

        def get_Xhat(Z, xr):
            W = self.get_W(Z, xr)
            return W.dot(self.U)

        if hasattr(self, 'Z'):
            Z = self.Z
        else:
            Z = np.zeros(self.z_shape)
        if not hasattr(self, 'Z_series'):
            self.Z_series = [np.zeros(self.z_shape)]
        self.get_Xhat = get_Xhat
        self.get_grad = self.grad_Z
        self.loss_series = self.Z_loss_series
        self.grad_series = self.Z_grad_series
        self.params_series = self.Z_series
        Z = self.nonlinear_proxGrad(Z, maxItr=maxItr, with_prints=with_prints, fa_update_freq=fa_update_freq,
                                    nonneg=nonneg, forceNorm=forceNorm)
        if with_prints:
            self.print_performance(Z, maxItr)
        self.Z = Z
        self.Z_loss_series = self.loss_series
        self.Z_grad_series = self.grad_series
        self.Z_series = self.params_series

    def update_U(self, maxItr=100, with_prints=False, fa_update_freq=10, nonneg=True, forceNorm=True):
        if with_prints:
            print('updating U')

        def get_Xhat(U, xr):
            W = self.get_W(self.Z, xr)
            return W.dot(U)

        self.get_Xhat = get_Xhat
        self.get_grad = self.grad_U
        self.loss_series = self.U_loss_series
        self.grad_series = self.U_grad_series
        if not hasattr(self, 'U_series'):
            self.U_series = [np.copy(self.U)]
        self.params_series = self.U_series
        U = self.nonlinear_proxGrad(self.U, maxItr=maxItr, with_prints=with_prints, fa_update_freq=fa_update_freq,
                                    nonneg=nonneg, forceNorm=forceNorm)
        if with_prints:
            self.print_performance(U, maxItr)
        self.U = U
        self.U_loss_series = self.loss_series
        self.U_grad_series = self.grad_series
        self.U_series = self.params_series

    def get_k_bounds(self):
        XU = []
        for xg in self.Xg:
            if len(xg) > 0:
                XU.append(xg.dot(self.U.T))
        XU = np.vstack(XU)
        return [(max(0.01, np.percentile(x, 0.5)), np.percentile(x, 99.5)) for x in XU.T]

    def update_fa(self, fmin_itrs=10):
        k, b, y0, x0 = self.k, self.b, self.y0, self.x0

        def partial_grad(fargs):
            Resid = []
            W_linear = []
            for xg, xr in zip(self.Xg[1:], self.Xr[:-1]):
                if len(xg) > 0:
                    w_linear, w, wgrad = self.get_W(self.Z, xr, func_args=fargs, both=True)
                    W_linear.append(w_linear)
                    xh = w.dot(self.U)
                    if not self.withCoupling:
                        Resid.append(ne.evaluate("xh - xg").dot(self.U.T) / xg.shape[0])
            return Resid, W_linear

        def get_loss(fargs):
            loss = 0
            for xg, xr in zip(self.Xg[1:], self.Xr[:-1]):
                if len(xg) > 0:
                    w = self.get_W(self.Z, xr, func_args=fargs)
                    xh = w.dot(self.U)
                    if not self.withCoupling:
                        loss += np.linalg.norm(ne.evaluate("xh - xg")) ** 2 / xg.shape[0]
            return loss

        def grad_k(k, b=b, y0=y0, x0=x0):
            farg = {'k': k, 'b': b, 'y0': y0, 'x0': x0}
            Resid, W_linear = partial_grad(farg)
            grad = np.zeros(len(k))
            for resid, w_linear in zip(Resid, W_linear):
                ex = ne.evaluate("exp(-(w_linear-x0)*b)")
                dk = ne.evaluate("(1 + ex*(k/y0-1))")
                bg = ne.evaluate("resid*(dk**-1 - (ex*k/y0)*dk**-2)")
                bg[np.isnan(bg)] = 0
                grad += bg.sum(0)
            return grad

        def func_k(k, b=b, y0=y0, x0=x0):
            farg = {'k': k, 'b': b, 'y0': y0, 'x0': x0}
            return get_loss(farg)

        if not np.isnan(grad_k(k).max()):
            bounds = [(0.01, 100) for _ in k]
            k, nfeval, rc = fmin_tnc(func_k, k, fprime=grad_k, messages=0, maxfun=fmin_itrs, bounds=bounds)

        def grad_y0(y0, k=k, b=b, x0=x0):
            farg = {'k': k, 'b': b, 'y0': y0, 'x0': x0}
            Resid, W_linear = partial_grad(farg)
            grad = np.zeros(len(y0))
            for resid, w_linear in zip(Resid, W_linear):
                ex = ne.evaluate("exp(-(w_linear-x0)*b)")
                dk = ne.evaluate("(1 + ex*(k/y0-1))")
                bg = ne.evaluate("resid*(ex*k/y0**2)*dk**-2")
                bg[np.isnan(bg)] = 0
                grad += bg.sum(0)
            return grad

        def func_y0(y0, k=k, b=b, x0=x0):
            farg = {'k': k, 'b': b, 'y0': y0, 'x0': x0}
            return get_loss(farg)

        if not np.isnan(grad_y0(y0).max()):
            bounds = [(1e-5, ki / 2.) for ki in k]
            y0, nfeval, rc = fmin_tnc(func_y0, y0, fprime=grad_y0, messages=0, maxfun=fmin_itrs, bounds=bounds)

        def grad_b(b, k=k, y0=y0, x0=x0):
            farg = {'k': k, 'b': b, 'y0': y0, 'x0': x0}
            Resid, W_linear = partial_grad(farg)
            grad = 0
            for resid, w_linear in zip(Resid, W_linear):
                ex = ne.evaluate("exp(-(w_linear-x0)*b)")
                dk = ne.evaluate("(1 + ex*(k/y0-1))")
                bg = ne.evaluate("resid*(ex*(w_linear-x0)*k*(k/y0-1))*dk**-2")
                bg[np.isnan(bg)] = 0
                grad += bg.sum()
            return grad

        def func_b(b, k=k, y0=y0, x0=x0):
            farg = {'k': k, 'b': b, 'y0': y0, 'x0': x0}
            return get_loss(farg)

        if not np.isnan(grad_b(b).max()):
            bounds = [(0.1, 6.)]
            b, nfeval, rc = fmin_tnc(func_b, b, fprime=grad_b, messages=0, maxfun=fmin_itrs, bounds=bounds)

        def grad_x0(x0, k=k, y0=y0, b=b):
            farg = {'k': k, 'b': b, 'y0': y0, 'x0': x0}
            Resid, W_linear = partial_grad(farg)
            grad = 0
            for resid, w_linear in zip(Resid, W_linear):
                ex = ne.evaluate("exp(-(w_linear-x0)*b)")
                dk = ne.evaluate("(1 + ex*(k/y0-1))")
                bg = ne.evaluate("resid*(ex*-b*k*(k/y0-1))*dk**-2")
                bg[np.isnan(bg)] = 0
                grad += bg.sum()
            return grad

        def func_x0(x0, k=k, y0=y0, b=b):
            farg = {'k': k, 'b': b, 'y0': y0, 'x0': x0}
            return get_loss(farg)

        if not np.isnan(grad_x0(x0).max()):
            bounds = [(-1, 1)]
            x0, nfeval, rc = fmin_tnc(func_x0, x0, fprime=grad_x0, messages=0, maxfun=fmin_itrs, bounds=bounds)
        self.set_fa(k, b, y0, x0)
        self.k, self.b, self.y0, self.x0 = k, b, y0, x0
