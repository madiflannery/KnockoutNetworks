import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;

/**
 * Created by Madison on 29/08/2015.
 */
public class ProcessFastDASOutput {

    public static final String NODE_FILE = "";
    public static final String NETWORK_FILE = "";
    public static final String FILE_PATH = "C:\\Users\\Madison\\Desktop\\FastDAS_RESULTS\\";
    public static final String[] FOLDERNAMES = {"outfile", "outfile1", "outfile2", "outfile3"};
    public static final String HICONF = "\\High_Conf\\";
    public static final int NODE_CYTO = 0;
    public static final int NODE_RDF = 2;



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

        for (String rdf : rdf_to_networks.keySet()) {
            ArrayList<String> files = rdf_to_networks.get(rdf);
            System.out.print(rdf + ": ");
            for(String file : files) {
                System.out.print(file + " || ");
            }
            System.out.print("\n");
        }
    }

    private static void processNetworkFile(HashMap<String, ArrayList<String>> rdf_to_networks, HashMap<String, String> cyto_to_rdf, String filename) {
        //open file
        File file = new File(filename);
        try(BufferedReader br = new BufferedReader(new FileReader(file))) {
            String line;

            //process the data
            while((line = br.readLine()) != null) {
                if(line.contains("id")) {
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
}
