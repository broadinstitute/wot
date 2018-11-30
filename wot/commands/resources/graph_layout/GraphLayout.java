import org.gephi.graph.api.Graph;
import org.gephi.graph.api.GraphController;
import org.gephi.graph.api.GraphModel;
import org.gephi.graph.api.Node;
import org.gephi.io.exporter.api.ExportController;
import org.gephi.io.importer.api.Container;
import org.gephi.io.importer.api.EdgeDirectionDefault;
import org.gephi.io.importer.api.ImportController;
import org.gephi.io.processor.plugin.DefaultProcessor;
import org.gephi.layout.plugin.forceAtlas2.ForceAtlas2;
import org.gephi.layout.plugin.fruchterman.FruchtermanReingold;
import org.gephi.layout.plugin.openord.OpenOrdLayout;
import org.gephi.layout.spi.Layout;
import org.gephi.project.api.ProjectController;
import org.gephi.project.api.Workspace;
import org.openide.util.Lookup;

import java.io.File;
import java.io.FileWriter;
import java.io.BufferedReader;
import java.io.FileReader;
import java.io.PrintWriter;
import java.io.IOException;
import java.util.Random;
import java.util.Set;
import java.util.List;
import java.util.Arrays;
import java.util.HashSet;
import java.util.HashMap;
import java.util.Map;


public class GraphLayout {

    private static void writeOutput(Graph g, boolean is3d, Set<String> formats, String output) {
        try {
            // ExporterCSV, ExporterDL, ExporterGDF, ExporterGEXF, ExporterGML, ExporterGraphML, ExporterPajek, ExporterVNA, PDFExporter, PNGExporter, SVGExporter
            ExportController ec = Lookup.getDefault().lookup(ExportController.class);
            for (String format : formats) {
                if (format.equals("txt")) {
                    PrintWriter pw = new PrintWriter(new FileWriter(output + (output.toLowerCase().endsWith("." + format) ? "" : "." + format)));
                    pw.print("id\tx\ty" + (is3d ? "\tz" : "") + "\n");
                    for (Node n : g.getNodes()) {
                        pw.print(n.getId());
                        pw.print("\t");
                        pw.print(n.X());
                        pw.print("\t");
                        pw.print(n.y());
                        if (is3d) {
                            pw.print("\t");
                            pw.print(n.z());
                        }
                        pw.print("\n");
                    }
                    pw.close();
                } else {
                    ec.exportFile(new File(output + (output.toLowerCase().endsWith("." + format) ? "" : "." + format)), ec.getExporter(format));
                }
            }
        } catch (IOException x) {
            x.printStackTrace();
            System.exit(1);
        }
    }

    public static void main(String[] args) throws IOException {
        File file = null;
        String output = null;
        String layoutString = "fa";
        int nsteps = 10000;
        Long seed = null;
        boolean barnesHutOptimize = true;
        int threadCount = Math.max(1, Runtime.getRuntime().availableProcessors() - 1);
        Double barnesHutTheta = null;
        Double jitterTolerance = null;
        Boolean linLogMode = null;
        Double scalingRatio = null;
        Boolean strongGravityMode = null;
        Double gravity = null;
        Boolean outboundAttractionDistribution = null;
        Integer barnesHutUpdateIter = null;
        Set<String> formats = new HashSet<>();
        File coordsFile = null;
        Boolean updateCenter = false;

        for (int i = 0; i < args.length; i++) {
            if (args[i].equals("--input")) {
                file = new File(args[++i]);
                if (!file.exists()) {
                    System.err.println(file + " not found.");
                    System.exit(1);
                }
            } else if (args[i].equals("--output")) {
                output = args[++i];
            } else if (args[i].equals("--layout")) {
                layoutString = args[++i];
            } else if (args[i].equals("--nsteps")) {
                nsteps = Integer.parseInt(args[++i]);
            } else if (args[i].equals("--barnesHutOptimize")) {
                barnesHutOptimize = args[++i].equalsIgnoreCase("true");
            } else if (args[i].equals("--threads")) {
                threadCount = Integer.parseInt(args[++i]);
            } else if (args[i].equals("--barnesHutTheta")) {
                barnesHutTheta = Double.parseDouble(args[++i]);
            } else if (args[i].equals("--jitterTolerance")) {
                jitterTolerance = Double.parseDouble(args[++i]);
            } else if (args[i].equals("--linLogMode")) {
                linLogMode = args[++i].equalsIgnoreCase("true");
            } else if (args[i].equals("--scalingRatio")) {
                scalingRatio = Double.parseDouble(args[++i]);
            } else if (args[i].equals("--gravity")) {
                gravity = Double.parseDouble(args[++i]);
            } else if (args[i].equals("--strongGravityMode")) {
                strongGravityMode = args[++i].equalsIgnoreCase("true");
            } else if (args[i].equals("--outboundAttractionDistribution")) {
                outboundAttractionDistribution = args[++i].equalsIgnoreCase("true");
            } else if (args[i].equals("--seed")) {
                seed = Long.parseLong(args[++i]);
            } else if (args[i].equals("--format")) {
                formats.add(args[++i]);
            } else if (args[i].equals("--barnesHutUpdateIter")) {
                barnesHutUpdateIter = Integer.parseInt(args[++i]);
            } else if (args[i].equals("--updateCenter")) {
                updateCenter = args[++i].equalsIgnoreCase("true");
            } else if (args[i].equals("--coords")) {
                coordsFile = new File(args[++i]);
                if (!coordsFile.exists()) {
                    System.err.println(coordsFile + " not found.");
                    System.exit(1);
                }
            } else {
                System.err.println("Unknown argument " + args[i]);
            }
        }
        if (formats.size() == 0) {
            formats.add("txt");
        }

        ProjectController pc = Lookup.getDefault().lookup(ProjectController.class);
        pc.newProject();
        Workspace workspace = pc.getCurrentWorkspace();
        ImportController importController = Lookup.getDefault().lookup(ImportController.class);
        GraphModel graphModel = Lookup.getDefault().lookup(GraphController.class).getGraphModel();
        Container container = importController.importFile(file);
        container.getLoader().setEdgeDefault(EdgeDirectionDefault.UNDIRECTED);
        Graph g = graphModel.getUndirectedGraph();
        Graph graph = graphModel.getGraph();
        importController.process(container, new DefaultProcessor(), workspace);

        Layout layout = null;
        if (layoutString.equals("fa")) {
            layout = new ForceAtlas2(null);
        } else if (layoutString.equals("fr")) {
            layout = new FruchtermanReingold(null);
        } else if (layoutString.equals("oo")) {
            layout = new OpenOrdLayout(null);
        } else if (layoutString.equals("fa_3d")) {
            layout = new org.gephi.layout.plugin.forceAtlas2_3d.ForceAtlas2(null);
        } else if (layoutString.equals("iso")) {
            layout = new org.alexandrebarao.isometriclayout.IsometricLayout(null);
        } else {
            System.err.println("Unknown layout " + layoutString);
            System.exit(1);
        }
        layout.setGraphModel(graphModel);

        Random random = seed != null ? new Random(seed) : new Random();
        final boolean is3d = layout instanceof org.gephi.layout.plugin.forceAtlas2_3d.ForceAtlas2;
        for (Node node : graph.getNodes()) {
            node.setX((float) ((0.01 + random.nextDouble()) * 1000) - 500);
            node.setY((float) ((0.01 + random.nextDouble()) * 1000) - 500);
            if (is3d) {
                node.setZ((float) ((0.01 + random.nextDouble()) * 1000) - 500);
            }
        }
        if (coordsFile != null) {
            Map<Object, Node> idToNode = new HashMap<>();
            for (Node n : g.getNodes()) {
                idToNode.put(n.getId(), n);

            }
            BufferedReader br = new BufferedReader(new FileReader(coordsFile));
            String sep = "\t";
            String s = br.readLine();
            for (String test : new String[]{"\t", ","}) {
                if (s.indexOf(test) != -1) {
                    sep = test;
                    break;
                }
            }
            List<String> header = Arrays.asList(s.split(sep));
            int idIndex = header.indexOf("id");
            int xIndex = header.indexOf("x");
            int yIndex = header.indexOf("y");
            int zIndex = header.indexOf("z");
            while ((s = br.readLine()) != null) {
                String[] tokens = s.split(sep);
                String id = tokens[idIndex];
                Node n = idToNode.get(id);
                if (n != null) {
                    n.setX(Float.parseFloat(tokens[xIndex]));
                    n.setY(Float.parseFloat(tokens[yIndex]));
                    if (is3d) {
                        n.setZ(Float.parseFloat(tokens[zIndex]));
                    }
                } else {
                    System.err.println(id + " not found");
                }
            }
            br.close();
        }

        if (layout instanceof ForceAtlas2) {
            ForceAtlas2 fa = (ForceAtlas2) layout;
            fa.setBarnesHutOptimize(barnesHutOptimize);
            if (barnesHutTheta != null) {
                fa.setBarnesHutTheta(barnesHutTheta);
            }
            if (jitterTolerance != null) {
                fa.setJitterTolerance(jitterTolerance);
            }
            if (linLogMode != null) {
                fa.setLinLogMode(linLogMode);
            }
            if (scalingRatio != null) {
                fa.setScalingRatio(scalingRatio);
            }
            if (strongGravityMode != null) {
                fa.setStrongGravityMode(strongGravityMode);
            }
            if (gravity != null) {
                fa.setGravity(gravity);
            }
            if (outboundAttractionDistribution != null) {
                fa.setOutboundAttractionDistribution(outboundAttractionDistribution);
            }
            fa.setThreadsCount(threadCount);
        }
        if (layout instanceof org.gephi.layout.plugin.forceAtlas2_3d.ForceAtlas2) {
            org.gephi.layout.plugin.forceAtlas2_3d.ForceAtlas2 fa = (org.gephi.layout.plugin.forceAtlas2_3d.ForceAtlas2) layout;
            fa.setBarnesHutOptimize(barnesHutOptimize);
            if (barnesHutTheta != null) {
                fa.setBarnesHutTheta(barnesHutTheta);
            }
            if (jitterTolerance != null) {
                fa.setJitterTolerance(jitterTolerance);
            }
            if (linLogMode != null) {
                fa.setLinLogMode(linLogMode);
            }
            if (scalingRatio != null) {
                fa.setScalingRatio(scalingRatio);
            }
            if (strongGravityMode != null) {
                fa.setStrongGravityMode(strongGravityMode);
            }
            if (gravity != null) {
                fa.setGravity(gravity);
            }
            if (outboundAttractionDistribution != null) {
                fa.setOutboundAttractionDistribution(outboundAttractionDistribution);
            }
            if (barnesHutUpdateIter != null) {
                fa.setUpdateBarnesHutIter(barnesHutUpdateIter);
            }
            if (updateCenter != null) {
                fa.setUpdateCenter(updateCenter);
            }
            fa.setThreadsCount(threadCount);
        }

        layout.initAlgo();
        final Set<String> _formats = formats;
        final String _output = output;
        final Graph _g = g;
        final Layout _layout = layout;
        Thread shutdownThread = new Thread() {
            @Override
            public void run() {
                _layout.endAlgo();
                writeOutput(_g, is3d, _formats, _output);
            }
        };
        Runtime.getRuntime().addShutdownHook(shutdownThread);
        int lastPercent = 0;
        for (int i = 0; i < nsteps; i++) {
            layout.goAlgo();
            int percent = (int) Math.floor(100 * (i + 1.0) / nsteps);
            if (percent != lastPercent) {
                System.out.print("*");
                lastPercent = percent;
                if (percent % 25 == 0) {
                    System.out.println(percent + "%");
                }
            }
        }
        if (nsteps > 0) {
            System.out.println();
        }
        layout.endAlgo();
//        Runtime.getRuntime().removeShutdownHook(shutdownThread);
//        writeOutput(g, is3d, formats, output);
    }
}
