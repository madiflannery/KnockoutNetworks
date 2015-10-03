import java.io.*;
import java.lang.reflect.Array;
import java.util.*;

/**
 * Created by Madison on 29/08/2015.
 */
public class ProcessFastDASOutput {

    public static final String NODE_FILE = "";
    public static final String FILE_PATH = "C:\\Users\\Madison\\Google Drive\\MSci\\Research Project\\KnockoutNetworks\\egf_old_results\\";
    public static final String OUTPUT_FILE = FILE_PATH + "output_processed_all_2.txt";
    public static final String OUTPUT_FILE_PHENO = FILE_PATH + "output_processed_all_2_pheno.txt";
    public static final String FASTDAS_COL_MAP_FILE = "C:\\Users\\Madison\\Google Drive\\MSci\\Research Project\\KnockoutNetworks\\egf_old_outfiles\\fastdas_col_mapping.txt";
    public static final String UNIPROT_MAPPING_FILE = "C:\\Users\\Madison\\Documents\\KnockoutNetworks\\input_files\\screen_gene_list_uniprot.txt";
    public static final String PHENOTYPE_DATA = "C:\\Users\\Madison\\Desktop\\Data Mapping 2\\PHENOTYPE_EGF.txt";
    public static final String[] FOLDERNAMES = {"outfile", "outfile1", "outfile2"}; //, "outfile3"};
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
        HashMap<String, String> gene_to_phenostr = read_phenotype_data_str(PHENOTYPE_DATA);
        HashMap<String, ArrayList<String>> gene_to_pheno = read_phenotype_data(PHENOTYPE_DATA);
        HashMap<String, String> uniprot_to_gene = read_uniprot_gene_mapping(UNIPROT_MAPPING_FILE, gene_to_phenostr);


        HashMap<String, HashSet<String>> gene_names_to_network_rdf = networkToGeneNames(rdf_to_networks, fastdas_col_mappings, uniprot_to_gene);
        ArrayList<String> rdf_list = new ArrayList<String>(rdf_to_networks.keySet());
        createOutputFile(gene_names_to_network_rdf, rdf_list, gene_to_pheno, false);
        createOutputFile(gene_names_to_network_rdf, rdf_list, gene_to_pheno, true);

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

    private static void createOutputFile(HashMap<String, HashSet<String>> gene_names_to_network_rdf, ArrayList<String> rdfs, HashMap<String, ArrayList<String>> gene_to_pheno, boolean phenotype_on) {
        ArrayList<ArrayList<String>> output = new ArrayList<>();
        ArrayList<String> ko_names = new ArrayList(gene_names_to_network_rdf.keySet());

        //build header line first
        ArrayList<String> temp_header = new ArrayList(gene_names_to_network_rdf.keySet());
        temp_header.add(0, "RDF");
        output.add(temp_header);

        if(phenotype_on) {
            //CELL COUNT
            ArrayList<String> temp_count = new ArrayList<>();
            temp_count.add("COUNT");
            ArrayList<String> temp_vim = new ArrayList<>();
            temp_vim.add("VIM");
            for (String ko : ko_names) {
                //for each KO
                //get the phenotype data
                //add to lists
                ArrayList<String> pheno = gene_to_pheno.get(ko.split(";")[0].toUpperCase().trim());
                if (pheno != null) {
                    temp_count.add(pheno.get(0));
                    temp_vim.add(pheno.get(1));
                } else {
                    temp_count.add("NA");
                    temp_vim.add("NA");
                }
            }
            output.add(temp_count);
            output.add(temp_vim);
        }

        //for each rdf
        for(String rdf : rdfs) {
            if(rdf != null) { // && rdf.toLowerCase().contains("reaction")) {
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

        File outfile;
        if(phenotype_on) {
            outfile = new File(OUTPUT_FILE_PHENO);
        }
        else {
            outfile = new File(OUTPUT_FILE);
        }
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
                String temp = file.replace(FILE_PATH, "");
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

    public static HashMap<String, String> read_uniprot_gene_mapping(String filename, HashMap<String, String> gene_to_pheno) {
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
                String uprot = s[1];
                String gene = s[0] + gene_to_pheno.get(s[0].toUpperCase().trim());
                uniprot_to_gene.put(uprot, gene);
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

    public static HashMap<String, String> read_phenotype_data_str(String filename) {
        //storage of data
        HashMap<String, String> gene_to_pheno = new HashMap<>();

        //open file
        File file = new File(filename);
        try(BufferedReader br = new BufferedReader(new FileReader(file))) {
            String line;
            line = br.readLine();
            //process the data
            while((line = br.readLine()) != null) {
                String[] s = line.split("\t");
                if(gene_to_pheno.get(s[0].trim()) != null){
                    System.out.println("DUPLICATE");
                }
                if(s.length == 1) {
                    gene_to_pheno.put(s[0].trim().toUpperCase(), " NODATA");
                }
                else {
                    if(Float.parseFloat(s[2]) >= 100) {
                        gene_to_pheno.put(s[0].trim().toUpperCase(), ";" + s[1] + ";UP");
                    }
                    else {
                        gene_to_pheno.put(s[0].trim().toUpperCase(), ";" + s[1] + ";DOWN");
                    }
                }
            }

            br.close();
        } catch(IOException e) {
            //exit on error
            System.err.println(filename);
            e.printStackTrace();
            System.exit(1);
        }
        return gene_to_pheno;
    }


    public static HashMap<String, ArrayList<String>> read_phenotype_data(String filename) {
        //storage of data
        HashMap<String, ArrayList<String>> gene_to_pheno = new HashMap<>();

        //open file
        File file = new File(filename);
        try(BufferedReader br = new BufferedReader(new FileReader(file))) {
            String line;
            line = br.readLine();
            //process the data
            while((line = br.readLine()) != null) {
                ArrayList<String> tmp = new ArrayList<>();
                String[] s = line.split("\t");
                if(gene_to_pheno.get(s[0].trim()) != null){
                    System.out.println("DUPLICATE");
                }
                if(s.length == 1) {
                    gene_to_pheno.put(s[0], null);
                }
                else {
                    tmp.add(s[1]);
                    tmp.add(s[2]);
                    gene_to_pheno.put(s[0].trim().toUpperCase(), tmp);
                }
            }

            br.close();
        } catch(IOException e) {
            //exit on error
            System.err.println(filename);
            e.printStackTrace();
            System.exit(1);
        }
        return gene_to_pheno;
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
