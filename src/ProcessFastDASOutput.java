import java.io.*;
import java.lang.reflect.Array;
import java.util.*;

/**
 * Created by Madison on 29/08/2015.
 */
public class ProcessFastDASOutput implements ProcessFastDAS_if{

    public static final String NODE_FILE = "";
    public static final String FILE_PATH = "C:\\Users\\Madison\\Google Drive\\MSci\\Research Project\\KnockoutNetworks\\normal_results\\";
    public static final String OUTPUT_FILE = FILE_PATH + "output_processed_onlyreacs.txt";
    public static final String FASTDAS_COL_MAP_FILE = "C:\\Users\\Madison\\Documents\\KnockoutNetworks\\fastdas_col_mapping.txt";
    public static final String UNIPROT_MAPPING_FILE = "C:\\Users\\Madison\\Documents\\KnockoutNetworks\\input_files\\screen_gene_list_uniprot.txt";
    public static final String[] FOLDERNAMES = {"outfile", "outfile1", "outfile2", "outfile3"};
    public static final String HICONF = "\\High_Conf\\";
    public static final int NODE_CYTO = 0;
    public static final int NODE_RDF = 2;
    public static int comp_count = 0;


    public static void main(String[] args) {

        HashMap<String, ArrayList<String>> rdf_to_networks = new HashMap<>();

        for(String folder_name : FOLDERNAMES) {
            //process node file
            HashMap<String, String> cyto_to_rdf = readNodeFile(FILE_PATH+folder_name+HICONF+findNodeFilename(folder_name));
            ArrayList<String> network_files = findNetworkFilenames(folder_name);
            //process network files
            for(String net_file : network_files) {
                processNetworkFile(rdf_to_networks, cyto_to_rdf, FILE_PATH+folder_name+HICONF+net_file);
            }
        }

        ArrayList<HashMap<String, String>> fastdas_col_mappings = read_fastdas_col_mapping(FASTDAS_COL_MAP_FILE);
        HashMap<String, String> uniprot_to_gene = read_uniprot_gene_mapping(UNIPROT_MAPPING_FILE);

        HashMap<String, HashSet<String>> gene_names_to_network_rdf = networkToGeneNames(rdf_to_networks, fastdas_col_mappings, uniprot_to_gene);
        ArrayList<String> rdf_list = new ArrayList<String>(rdf_to_networks.keySet());
        createOutputFile(gene_names_to_network_rdf, rdf_list);

        int numgreater =0;
        for (String rdf : rdf_to_networks.keySet()) {
            ArrayList<String> files = rdf_to_networks.get(rdf);
//            System.out.print(rdf + ": ");
//            System.out.print(files.size());
            if(files.size()>1) {numgreater++;}
//            for(String file : files) {
//                System.out.print(file + " || ");
//            }
//            System.out.print("\n");
        }
        System.out.println("\nNUM NON-EMPTY COMPARISONS: " + comp_count);
        System.out.println("\nNUM RDFS: " + rdf_to_networks.keySet().size());
        System.out.println("\nNUM RDFS GREATER THAN ONE: " + numgreater);
    }

    private static void createOutputFile(HashMap<String, HashSet<String>> gene_names_to_network_rdf, ArrayList<String> rdfs) {
        ArrayList<ArrayList<String>> output = new ArrayList<>();
        ArrayList<String> ko_names = new ArrayList(gene_names_to_network_rdf.keySet());

        //build header line first
        ArrayList<String> temp_header = new ArrayList(gene_names_to_network_rdf.keySet());
        temp_header.add(0, "RDF");
        output.add(temp_header);

        //for each rdf
        for(String rdf : rdfs) {
            if(rdf.toLowerCase().contains("reaction")) {
                ArrayList<String> temp_line = new ArrayList<>();
                temp_line.add(rdf);
                //for each ko
                for (String ko : ko_names) {
                    //if rdf in KO's map
                    if (gene_names_to_network_rdf.get(ko).contains(rdf)) {
                        //put 1 in output
                        temp_line.add("1");
                    } else {//else
                        //put 0 in output
                        temp_line.add("0");
                    }
                }
                output.add(temp_line);
            }
        }

        File outfile = new File(OUTPUT_FILE);
        try {
            //open normal file and active file
            PrintWriter out = new PrintWriter(new BufferedWriter(new FileWriter(outfile)));

            for (ArrayList<String> line : output) {
                StringJoiner strbuild = new StringJoiner("\t");
                for (String s : line) {
                    strbuild.add(s);
                }
                out.println(strbuild.toString());
            }
            out.close();
        }
        catch (Exception e) {
            System.err.println(e.toString());
        }

    }

    public static HashMap<String, HashSet<String>> networkToGeneNames(HashMap<String, ArrayList<String>> rdf_to_networks, ArrayList<HashMap<String, String>> fastdas_col_mappings, HashMap<String, String> uniprot_to_gene) {
        HashMap<String, HashSet<String>> gene_names_to_network_rdf = new HashMap<>();

        for (String rdf : rdf_to_networks.keySet()) {
            ArrayList<String> files = rdf_to_networks.get(rdf);
            for (String file : files) {
                //get name of folder and the number of comparison
                String temp = file.replace("C:\\Users\\Madison\\Google Drive\\MSci\\Research Project\\KnockoutNetworks\\normal_results\\", "");
                temp = temp.replace("High_Conf\\Network_Control-", "");
                temp = temp.substring(0, temp.indexOf("_"));
                //System.out.print(temp);
                int numfolder = Arrays.asList(FOLDERNAMES).indexOf(temp.split("\\\\")[0].trim());
                String numcomp = temp.split("\\\\")[1].trim();
                String uniprot = fastdas_col_mappings.get(numfolder).get(numcomp);
                String genename = uniprot_to_gene.get(uniprot);
                //System.out.print(" - " + uniprot + ", " + genename + "\n");
                if(genename == null) {
                    System.out.println("ERROR: UNIPROT DID NOT GIVE GENE NAME");
                }
                if(gene_names_to_network_rdf.get(genename) == null) {
                    //add to hash map
                    HashSet<String> hs = new HashSet<>();
                    hs.add(rdf);
                    gene_names_to_network_rdf.put(genename, hs);
                }
                else {
                    gene_names_to_network_rdf.get(genename).add(rdf);
                }
            }
        }
        return gene_names_to_network_rdf;

    }

    public static HashMap<String, String> read_uniprot_gene_mapping(String filename) {
        //storage of data
        HashMap<String, String> uniprot_to_gene = new HashMap<>();

        //open file
        File file = new File(filename);
        try(BufferedReader br = new BufferedReader(new FileReader(file))) {
            String line;

            //process the data
            while((line = br.readLine()) != null) {
                String[] s = line.split("\t");
                if(uniprot_to_gene.get(s[1]) != null){
                    System.out.println("DUPLICATE");
                }
                uniprot_to_gene.put(s[1], s[0]);
            }

            br.close();
        } catch(IOException e) {
            //exit on error
            System.err.println(filename);
            e.printStackTrace();
            System.exit(1);
        }
        return uniprot_to_gene;
    }

    public static void processNetworkFile(HashMap<String, ArrayList<String>> rdf_to_networks, HashMap<String, String> cyto_to_rdf, String filename) {
        //open file
        File file = new File(filename);
        try(BufferedReader br = new BufferedReader(new FileReader(file))) {
            String line;
            boolean found = false;
            //process the data
            while((line = br.readLine()) != null) {
                if(line.contains("id")) {
                    found = true;
                    String cyto_num = line.split(" ")[line.split(" ").length-1];
                    if(rdf_to_networks.get(cyto_to_rdf.get(cyto_num)) == null) {
                        //DOESNT EXIST
                        ArrayList<String> toadd = new ArrayList<>();
                        toadd.add(filename);
                        rdf_to_networks.put(cyto_to_rdf.get(cyto_num), toadd);
                    }
                    else { //EXISTS
                        rdf_to_networks.get(cyto_to_rdf.get(cyto_num)).add(filename);
                    }
                }
            }

            br.close();

            if(found) {
                comp_count++;
                //System.out.println(filename);
            }
        } catch(IOException e) {
            //exit on error
            System.err.println(filename);
            e.printStackTrace();
            System.exit(1);
        }
    }

    public static String findNodeFilename(String folder_name) {
        File folder = new File(FILE_PATH+folder_name+HICONF);
        File[] listOfFiles = folder.listFiles();

        for (int i = 0; i < listOfFiles.length; i++) {
            if (listOfFiles[i].isFile() && listOfFiles[i].getName().contains("NodeAttributes")) {
                return (listOfFiles[i].getName());
            }
        }
        return null;
    }

    public static ArrayList<String> findNetworkFilenames(String folder_name) {
        ArrayList<String> res = new ArrayList<>();
        File folder = new File(FILE_PATH+folder_name+HICONF);
        File[] listOfFiles = folder.listFiles();

        for (int i = 0; i < listOfFiles.length; i++) {
            if (listOfFiles[i].isFile() && listOfFiles[i].getName().contains("gml")) {
                res.add(listOfFiles[i].getName());
            }
        }
        return res;
    }


    public static HashMap<String, String> readNodeFile(String filename) {
        //storage of data
        HashMap<String, String> cyto_to_rdf = new HashMap<>();

        //open file
        File file = new File(filename);
        try(BufferedReader br = new BufferedReader(new FileReader(file))) {
            String line;

            //process the data
            while((line = br.readLine()) != null) {
                String[] s = line.split("\t");
                if(cyto_to_rdf.get(s[NODE_CYTO]) != null){
                    System.out.println("DUPLICATE");
                }
                cyto_to_rdf.put(s[NODE_CYTO], s[NODE_RDF]);
            }

            br.close();
        } catch(IOException e) {
            //exit on error
            System.err.println(filename);
            e.printStackTrace();
            System.exit(1);
        }
        return cyto_to_rdf;
    }


    public static ArrayList<HashMap<String, String>> read_fastdas_col_mapping(String filename) {
        ArrayList<HashMap<String, String>> fastdas_col_mappings = new ArrayList<>();
        //open file
        File file = new File(filename);
        try(BufferedReader br = new BufferedReader(new FileReader(file))) {
            String line;
            HashMap<String, String> hm = null;
            //process the data
            while((line = br.readLine()) != null) {
                if(line.contains("----")) {
                    if(hm!= null) {
                        fastdas_col_mappings.add(new HashMap<String, String>(hm));
                    }
                    //new hash map
                    hm = new HashMap<>();
                    //TODO: WILL THIS WORK?
                }
                else {
                    String[] s = line.split(", ");

                    hm.put(s[0], s[1]);
                }
            }
            fastdas_col_mappings.add(new HashMap<String, String>(hm));
            br.close();
        } catch(IOException e) {
            //exit on error
            System.err.println(filename);
            e.printStackTrace();
            System.exit(1);
        }

        return fastdas_col_mappings;
    }
}
