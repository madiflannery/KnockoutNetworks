import java.lang.reflect.Array;
import java.util.regex.Pattern;
import org.biopax.paxtools.model.level3.*;
import org.biopax.paxtools.model.Model;

import java.io.*;
import java.util.*;
import java.util.regex.*;

import org.biopax.paxtools.impl.level3.ModificationFeatureImpl;
import org.biopax.paxtools.io.SimpleIOHandler;

/**
 * Created by Madison on 8/08/2015.
 */
public class KnockoutNetworks {
    //The biopax model
    private static Model biopaxModel = null;

    //UniProt Accession Patterns:
    //Pattern 1 : [A-N,R-Z][0-9][A-Z][A-Z, 0-9][A-Z, 0-9][0-9]
    private static final Pattern UNIPROT_PATTERN_ONE = Pattern.compile("([A-N,R-Z][0-9][A-Z][A-Z,0-9][A-Z, 0-9][0-9])");
    //Pattern 2 : [O,P,Q][0-9][A-Z, 0-9][A-Z, 0-9][A-Z, 0-9][0-9]
    private static final Pattern UNIPROT_PATTERN_TWO = Pattern.compile("([O,P,Q][0-9][A-Z,0-9][A-Z,0-9][A-Z," +
            "0-9][0-9])");
    private static final Pattern UNIPROT_PATTERN_THREE = Pattern.compile("[O,P,Q][0-9][A-Z,0-9]{3}[0-9]|[A-N,R-Z][0-9]([A-Z][A-Z,0-9]{2}[0-9]){1,2}");

    //This is the comment flag used to denote a species is a set.
    private static final String REACTOME_ENTITYSET_FLAG = "Converted from EntitySet in Reactome";

    //Filenames
    private static final String PARAMS_FILE = "parameters.txt";
    private static final String NETWORK_FILE = "species.txt";
    private static final String KNOCKOUT_GENE_LIST = "screen_gene_list_uniprot.txt";
    private static final String DIFF_EXP_FILE = "off_genes_fromdiff_uniprot.txt";
    private static final String OUTPUT_FILENAME = "outfile";
    private static String biopaxFilename;

    //Directories
    private static String outputDir;
    private static String inputDir;

    private static final int RDF_COL = 1; //0 indexed, for network species.txt file

    private static int rdfLen = 0;

    /**
     * Main function
     * @param args Command line arguments
     */
    public static void main(String[] args) {
        //read in params file & process
        System.out.println("Initialising parameters...");
        init();

        //open and process all files
        System.out.println("Reading data...");
        ArrayList<ArrayList<String>> network = readFile(inputDir + NETWORK_FILE);
        File biopaxFile = new File(inputDir + biopaxFilename);


        System.out.println("Deleting old output files...");
        // Delete all output files if they exist
        deleteAllFiles();

        System.out.println("Creating UniProt Map...");
        HashMap<String, HashSet<Protein>> uniprotMap = getUniProtMap(biopaxFile);

        //get the length of the garbage at the start of the RDF's
        rdfLen = (uniprotMap.values().iterator().next().iterator().next().getRDFId().indexOf("#"))+1;

        //BACKGROUND FILE PROCESSING
        ProcessedData background = new ProcessedData();
        background.setGeneList(openGeneList(inputDir + DIFF_EXP_FILE));
        background.setProteins(processBackgroundGeneList(uniprotMap, background.getGeneList()));
        background.setComplexes(mapAllComplexes(background.getProteins()));

        //FUNCTIONAL SCREEN PROCESSING
        ArrayList<String> geneList_screen = openGeneList(inputDir + KNOCKOUT_GENE_LIST);
        ArrayList<ProcessedData> knockouts = new ArrayList<>();
        for(String uniprot : geneList_screen) {
            ProcessedData ko = new ProcessedData();
            ArrayList<String> uprot = new ArrayList<>();
            uprot.add(uniprot);
            ko.setGeneList(uprot);
            ko.setProteins(processBackgroundGeneList(uniprotMap, ko.getGeneList()));
            ko.setComplexes(mapAllComplexes(ko.getProteins()));
            knockouts.add(ko);
        }

        //CREATE FILES
        HashMap<String, String> bg_network_rdf = createBackgroundRdfs(network, background.getProteins(), background.getComplexes());
        addBGtoOutputFile(network, bg_network_rdf);
        determineToMap(network, bg_network_rdf, background, knockouts);
        addKOtoOutputFile(network, bg_network_rdf, background, knockouts);
        //System.out.print(network);
    }



    /**
     * Initialise the necessary parameters using the parameters.txt file
     */
    private static void init() {
        //temp storage of data
        HashMap<String, String> params = new HashMap<>();

        //read all lines from file
        File file = new File(PARAMS_FILE);
        try (BufferedReader br = new BufferedReader(new FileReader(file))) {
            String line;
            //split on '=' to separate variable and value
            while ((line = br.readLine()) != null) {
                String[] splitted = line.split("=");
                params.put(splitted[0].trim(), splitted[1].trim());
            }
            br.close();
        } catch (IOException e) {
            //exit on error
            e.printStackTrace();
            System.exit(1);
        }

        //set the parameters
        biopaxFilename = params.get("biopax_file");

        //if theres no \ at the end of the filepath
        if(params.get("out_file_directory").substring(params.get("out_file_directory").length()).equals("\\")) {
            outputDir = params.get("out_file_directory");
        }
        else {
            outputDir = params.get("out_file_directory") + "\\";
        }

        //if theres no \ at the end of the filepath
        if(params.get("in_file_directory").substring(params.get("in_file_directory").length()).equals("\\")) {
            inputDir = params.get("in_file_directory");
        }
        else {
            inputDir = params.get("in_file_directory") + "\\";
        }

    }

    /**
     * Read a file.
     *
     * @param filename Filename of the file to read
     * @return ArrayList containing each line as ArrayList of strings
     */
    public static ArrayList<ArrayList<String>> readFile(String filename) {
        //storage of data
        ArrayList<ArrayList<String>> finalData = new ArrayList<>();

        //open file
        File file = new File(filename);
        try(BufferedReader br = new BufferedReader(new FileReader(file))) {
            String line;

            //process the data
            while((line = br.readLine()) != null) {
                finalData.add(arrayToArrayList(line.split("\t")));
            }

            br.close();
        } catch(IOException e) {
            //exit on error
            System.err.println(filename);
            e.printStackTrace();
            System.exit(1);
        }

        return finalData;
    }

    private static ArrayList<String> arrayToArrayList(String[] arr) {
        ArrayList<String> arrlist = new ArrayList<>();
        for(String s : arr) {
            arrlist.add(s.trim());
        }
        return arrlist;
    }

    /**
     * Deletes all output files if they exist to avoid file corruption
     */
    private static void deleteAllFiles() {
        String temp_filename =OUTPUT_FILENAME + ".txt";
        deleteFile(temp_filename);
        int temp = 1;

        //Keep trying to create directory with new number until successful
        do{
            temp_filename = OUTPUT_FILENAME + "(" + temp + ")" + ".txt";
            temp++;
        } while(deleteFile(temp_filename));

    }

    /**
     * Deletes an individual file if it exists
     * @param filename The name of the file to delete
     * @return true if deleted, false if doesnt exist
     */
    private static boolean deleteFile(String filename) {
        //check file exists
        File file = new File(filename);
        boolean exists = file.exists();

        //delete file
        if(exists) {
            file.delete();
            return true;
        }
        else {
            return false;
        }
    }

    /**
     * Method for getting a UniProt map from the BioPAX model file.
     * Written by Liam
     * @param biopaxFile The reactome .owl file (as a File object)
     * @return HashMap where key = Uniprot ID, value = corresponding biopax Protein object(s)
     */
    private static HashMap<String, HashSet<Protein>> getUniProtMap(File biopaxFile) {
        try {
            //Read BioPAX file.
            biopaxModel = readModelFromFile(biopaxFile);
            //Map for return....
            HashMap<String, HashSet<Protein>> uniProtMap = new HashMap<>();
            //For each Protein in the model:
            for (Protein protein : biopaxModel.getObjects(Protein.class)) {
                //Get the UniProt accession.
                String uniProt = getUniProt(protein);
                //Warn if something went bad during getting UniProt accession.
                if (uniProt == null || uniProt.isEmpty()) {
                    if (protein.getComment().contains(REACTOME_ENTITYSET_FLAG)) {
                        StringBuilder warnString = new StringBuilder();
                        warnString.append("EntitySet skipped : ");
                        warnString.append(protein.getDisplayName());
                        warnString.append("( ");
                        warnString.append(protein.getRDFId().replaceAll(biopaxModel.getXmlBase(), ""));
                        warnString.append(")");
                        System.err.println(warnString.toString());
                        continue;
                    } else {
                        System.err.println("ERROR : " + protein.getRDFId().replaceAll(biopaxModel.getXmlBase(), ""));
                        continue;
                    }
                    //Otherwise, we got the Uniprot OK, check to see if we need to make a new entry in the map...
                } else if (uniProtMap.containsKey(uniProt)) {
                    //Have an entry, just append.
                    uniProtMap.get(uniProt).add(protein);
                } else {
                    //New entry, instantiate new set.
                    HashSet<Protein> newSet = new HashSet<>();
                    newSet.add(protein);
                    uniProtMap.put(uniProt, newSet);
                }
            }
            //Return uniprot map.
            return uniProtMap;
            //Pokemon handling - if -anything- throws an exception we want to spit the dummy.
        } catch (Exception e) {
            e.printStackTrace();
            System.exit(-1);
        }
        return null;
    }

    /**
     * Reads a BioPAX Level 3 file into memory, returning a paxtools Model object.
     * Written by Liam
     * @param biopaxFile a java.io.File instance containing BioPAX Level 3 data (e.g. new File("homo sapiens.owl"))
     * @return paxtools Model object
     */
    @SuppressWarnings("finally")
    private static Model readModelFromFile(File biopaxFile) {
        FileInputStream fis = null;
        Model model = null;
        try {
            SimpleIOHandler ioHandler = new SimpleIOHandler();
            fis = new FileInputStream(biopaxFile);
            model = ioHandler.convertFromOWL(fis);
        } catch (Exception e) {
            if (biopaxFile == null) {
                System.err.println("Read of BioPAX model failed - null file, check your file selection/file contents.");
                System.err.println("Shutting down.");
                System.exit(-1);
            } else {
                String filename = biopaxFile.getName();
                System.err.println("Read of BioPAX model failed - check " + filename);
                System.err.println("Shutting down.");
                e.printStackTrace();
                System.exit(-1);
            }
        } finally {
            if (fis != null) {
                try {
                    fis.close();
                } catch (IOException e) {
                    System.err.println("Shutting down.");
                    e.printStackTrace();
                    System.exit(-1);
                }
            }
            return model;
        }
    }

    /**
     * Gets the Uniprot ID associated with a Protein object. Note that there
     * are two Uniprot namespaces (see class variables). Most of the time,
     * uniprot IDs are stored in the Protein's 'entityReference' object.
     * If not, they're stored in the Xrefs for the object - this checks the
     * entityReference, then the Xrefs.
     * Written by Liam
     * @param protein a paxtools Protein object.
     * @return String containing the UniProt ID
     */
    private static String getUniProt(Protein protein) {
        String uniProt = "";
        EntityReference entityReference = protein.getEntityReference();
        if (entityReference != null) {
            //Check for uniprot in EntityReference.
            for (String s : entityReference.getName()) {
                Matcher matcher = UNIPROT_PATTERN_ONE.matcher(s);
                if (matcher.find()) {
                    uniProt = matcher.group(1);
                    break;
                } else {
                    matcher = UNIPROT_PATTERN_TWO.matcher(s);
                    if (matcher.find()) {
                        uniProt = matcher.group(1);
                        break;
                    }
                    else {
                        matcher = UNIPROT_PATTERN_THREE.matcher(s);
                        if (matcher.find()) {
                            uniProt = matcher.group(1);
                            break;
                        }
                    }
                }
            }
            //Check the Xrefs if the EntityReference didn't get a Uniprot.
            if (uniProt.equals("")) {
                if(protein.getXref() != null) {
                    HashSet<Xref> xrefs = new HashSet<Xref>(protein.getXref());
                    if(!xrefs.isEmpty()) {
                        //if (xrefs != null) {
                        for (Xref xref : xrefs) {
                            if (xref.getDb().equals("UniProt")) {
                                uniProt = xref.getId();
                                break;
                            }
                        }
                        System.err.println(protein.getRDFId() + " has no UniProts");
                    }
                }
                else {
                    System.err.println(protein.getRDFId() + " has no Xrefs");
                }
            }

        } else {
            System.err.println(protein.getRDFId() + " has no EntityReference");
        }


        return uniProt;
    }

    /**
     * Opens a file containing a list of genes and associated Uniprots of
     * the form gene \t uniprot
     * returns an ArrayList containing the list of uniprot ID's
     * @param filename the filename of the gene list
     * @return ArrayList<String> containing the list of uniprot ID's
     */
    private static ArrayList<String> openGeneList(String filename) {
        //storage of data
        ArrayList<String> finalData = new ArrayList<>();

        //open file
        File file = new File(filename);
        try(BufferedReader br = new BufferedReader(new FileReader(file))) {
            String line;

            //process the data
            while((line = br.readLine()) != null) {
                //System.out.println(line.split("\t")[1]);
                finalData.add(line.split("\t")[1]);
            }
            br.close();
        } catch(IOException e) {
            //exit on error
            e.printStackTrace();
            System.exit(1);
        }
        return finalData;
    }

    /**
     * Takes a list of uniprot ID's. for each uniprot ID, find all associated RDFs and
     * create a RDF -> protein mapping.
     * @param uniprotMap the reactome mapping from uniprot to proteins
     * @param geneList ArrayList<String> list of uniprot ID's
     * @return HashMap<String, Protein> containing a mapping from <RDF, PROTEIN>
     */
    private static HashMap<String, Protein> processBackgroundGeneList(HashMap<String, HashSet<Protein>> uniprotMap, ArrayList<String> geneList) {
        HashMap<String, Protein> mapped = new HashMap<>();
        //int numUniprotsMapped = 0;
        //int numProteinObjs = 0;

        //for each uniprot in the gene list
        for (String uniprot : geneList) {
            //get the protein objects from hashmap- upper case
            HashSet<Protein> temp = uniprotMap.get(uniprot.toUpperCase());
            //if we have a match
            if(temp != null) {
                //numUniprotsMapped++;

                //for each corresponding protein, add to map
                for(Protein p : temp) {
                    mapped.put(p.getRDFId().substring(rdfLen), p);
                    //numProteinObjs++;
                }
            }


        }
        //System.out.println("NUM UNIPROTS MAPPED: " + numUniprotsMapped);
        //System.out.println("NUM PROTEIN OBJECTS: " + numProteinObjs);
        return mapped;
    }

    /**
     * Creates a 2d array containing all traversals of complex trees for each protein.
     * @param rdf_proteins HashMap containing a mapping from RDF -> protein
     * @return ArrayList<ArrayList<PhysicalEntity>> list of all tree traversals.
     */
    private static ArrayList<ArrayList<PhysicalEntity>> mapAllComplexes(HashMap<String, Protein> rdf_proteins) {
        ArrayList<ArrayList<PhysicalEntity>> complexes = new ArrayList<>();

        //for each mapped phosphoprotein
        for(Protein p : rdf_proteins.values()) {
            //find the complex(es) it is associated with
            ArrayList<PhysicalEntity> comp = mapComplex(p);
            if(comp != null) { //if we found something, add to list
                complexes.add(comp);
            }
        }

        return complexes;
    }

    /**
     * Builds a list by traversing the complex tree from the leaf node (protein)
     * @param p A paxtools Protein object
     * @return A list of PhysicalEntity nodes
     */
    private static ArrayList<PhysicalEntity> mapComplex(Protein p) {
        //setup
        ArrayList<PhysicalEntity> seen = new ArrayList<>();
        LinkedList<PhysicalEntity> q = new LinkedList<>(); //queue
        q.add(p);
        seen.add(p);

        //if protein in no complexes
        if(p.getComponentOf().size() == 0) {
            return null;
        }

        //process the whole subtree above where we start
        while(!q.isEmpty()) {
            PhysicalEntity current = q.poll();
            for(PhysicalEntity ent : current.getComponentOf()) {
                if(!seen.contains(ent)) {
                    q.add(ent);
                    seen.add(ent);
                }
            }
        }
        return seen;
    }

    private static HashMap<String, String> createBackgroundRdfs(ArrayList<ArrayList<String>> network, HashMap<String, Protein> bg_prots, ArrayList<ArrayList<PhysicalEntity>> bg_complexes) {
        HashMap<String, String> rdf_to_status = new HashMap<>();
        for(ArrayList<String> line : network) {
            boolean setFree = true;
            if(line.get(RDF_COL).contains("Protein")) {
                Protein p = bg_prots.get(line.get(RDF_COL));
                if(p != null) {
                    //we have a matching protein, set inactive
                    rdf_to_status.put(line.get(RDF_COL), "INACTIVE");
                    setFree = false;
                }
            }
            else if(line.get(RDF_COL).contains("Complex")) {
                ArrayList<PhysicalEntity> complexTree = findComplexTree(line.get(RDF_COL), bg_complexes);
                if(complexTree != null) {
                    //we have a complex in tree, set inactive
                    rdf_to_status.put(line.get(RDF_COL), "INACTIVE");
                    setFree = false;
                }
            }
            if(setFree) {
                rdf_to_status.put(line.get(RDF_COL), "FREE");
            }
        }
        return rdf_to_status;
    }

    private static void addBGtoOutputFile(ArrayList<ArrayList<String>> network, HashMap<String, String> bg_network_rdf) {
        for(int i = 0; i < network.size(); i++) {
            network.get(i).add(bg_network_rdf.get(network.get(i).get(RDF_COL)));
        }
    }

    private static void addKOtoOutputFile(ArrayList<ArrayList<String>> network, HashMap<String, String> bg_network_rdf, ProcessedData background, ArrayList<ProcessedData> knockouts) {
        int count = 1;
        int nummapped = 0;
        for(ProcessedData ko : knockouts) {
            if(ko.getToAdd()) {
                boolean mapped = false;
                for (int i = 0; i < network.size(); i++) {
                    boolean setFree = true;
                    if (network.get(i).get(RDF_COL).contains("Protein")) {
                        Protein p = ko.getProteins().get(network.get(i).get(RDF_COL));
                        if (p != null) {
                            //we have a matching protein, set inactive
                            network.get(i).add("INACTIVE");
                            setFree = false;
                            mapped = true;
                        }
                    } else if (network.get(i).get(RDF_COL).contains("Complex")) {
                        ArrayList<PhysicalEntity> complexTree = findComplexTree(network.get(i).get(RDF_COL), ko.getComplexes());
                        if (complexTree != null) {
                            //we have a complex in tree, set inactive
                            network.get(i).add("INACTIVE");
                            setFree = false;
                            mapped = true;
                        }
                    }
                    if (setFree) {
                        network.get(i).add(bg_network_rdf.get(network.get(i).get(RDF_COL)));
                    }
                }
                if (mapped) {
                    nummapped++;
                }
                count++;
                if ((count % 50) == 0) {
                    //write current network list to file
                    writeToFile(network);

                    //reset network list
                    network = readFile(inputDir + NETWORK_FILE);
                    count++;

                    //re-add background RDFs
                    addBGtoOutputFile(network, bg_network_rdf);
                }
            }
        }
        //write current network list to file
        writeToFile(network);
        System.out.println(nummapped);
    }


    private static void determineToMap(ArrayList<ArrayList<String>> network, HashMap<String, String> bg_network_rdf, ProcessedData background, ArrayList<ProcessedData> knockouts) {
        for(ProcessedData ko : knockouts) {
            for(int i = 0; i < network.size(); i++) {
                boolean setFree = true;
                if(network.get(i).get(RDF_COL).contains("Protein")) {
                    Protein p = ko.getProteins().get(network.get(i).get(RDF_COL));
                    if(p != null) {
                        //we have a matching protein, set inactive
                        ko.setToAdd(true);
                    }
                }
                else if(network.get(i).get(RDF_COL).contains("Complex")) {
                    ArrayList<PhysicalEntity> complexTree = findComplexTree(network.get(i).get(RDF_COL), ko.getComplexes());
                    if(complexTree != null) {
                        //we have a complex in tree, set inactive
                        ko.setToAdd(true);
                    }
                }
            }
        }
    }

    private static void writeToFile(ArrayList<ArrayList<String>> network) {
        String filename = findFilename();
        File outfile = new File(filename);
        try {
            //open normal file and active file
            PrintWriter out = new PrintWriter(new BufferedWriter(new FileWriter(outfile)));

            for (ArrayList<String> line : network) {
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

    /**
     * Function to find a complex tree that a given complex is a part of. Returns null if not found.
     * @param rdf The RDF of the complex to search for
     * @param complexes ArrayList of complex trees to search
     * @return ArrayList - the complex tree containing complex. Null if no trees found
     */
    private static ArrayList<PhysicalEntity> findComplexTree(String rdf, ArrayList<ArrayList<PhysicalEntity>> complexes) {
        //for each tree
        for (ArrayList<PhysicalEntity> tree : complexes) {
            //for each physical entity in complex
            for(PhysicalEntity pe : tree) {
                //if we have a match, return tree
                if(pe.getRDFId().substring(rdfLen).equals(rdf)) {
                    return tree;
                }
            }
        }
        //not found
        return null;
    }

    private static String findFilename() {
        String temp_filename =OUTPUT_FILENAME + ".txt";
        File temp_file = new File(temp_filename);

        //If the directory doesnt exist, create it
        if(temp_file.exists()) {
            int temp = 1;

            //Keep trying to create directory with new number until successful
            do{
                temp_filename = OUTPUT_FILENAME + "(" + temp + ")" + ".txt";
                temp_file = new File(temp_filename);
                temp++;
            } while(temp_file.exists());
        }
        return temp_filename;
    }
}

