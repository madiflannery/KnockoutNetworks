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
    private static String biopaxFilename;

    //Directories
    private static final String IFILE_DIR = "input_files\\";
    private static String outputDir;

    private static final int RDF_COL = 1; //0 indexed
    private static final int INPUT_COL = 5;

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
        ArrayList<String[]> network = readFile(NETWORK_FILE, false);
        File biopaxFile = new File(biopaxFilename);


        System.out.println("Deleting old output files...");
        // Delete all output files if they exist
        deleteAllFiles();

        System.out.println("Creating UniProt Map...");
        HashMap<String, HashSet<Protein>> uniprotMap = getUniProtMap(biopaxFile);

        //get the length of the garbage at the start of the RDF's
        rdfLen = (uniprotMap.values().iterator().next().iterator().next().getRDFId().indexOf("#"))+1;

        //BACKGROUND FILE PROCESSING
        ProcessedData background = new ProcessedData();
        background.setGeneList(openGeneList(DIFF_EXP_FILE));
        background.setProteins(processBackgroundGeneList(uniprotMap, background.getGeneList()));
        background.setComplexes(mapAllComplexes(background.getProteins()));

        HashMap<String, String> rdfToBGVal = genBackgroundFile(network, bg_prots, bg_complexes);

        //FUNCTIONAL SCREEN PROCESSING
        ArrayList<String> geneList_screen = openGeneList(KNOCKOUT_GENE_LIST);
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
        if(params.get("file_directory").substring(params.get("file_directory").length()).equals("\\")) {
            outputDir = params.get("file_directory");
        }
        else {
            outputDir = params.get("file_directory") + "\\";
        }

    }

    //TODO: CHANGE THIS FUNCTION TO RETURN AN ARRAYLIST<ARRAYLIST<STRING>>
    //TODO: CHANGE THIS FUNCTION TO REMOVE ISPHOSPHO PARAMETER
    /**
     * Read the network or phosphoproteomics file.
     * If phosphoproteomics file, the file must have a header to process
     * to get cell line names.
     *
     * @param filename Filename of the file to read
     * @param isPhospho <code>True</code> if we need to process the cell line names
     * @return ArrayList containing each line as array of strings
     */
    public static ArrayList<String[]> readFile(String filename, boolean isPhospho) {
        //storage of data
        ArrayList<String[]> finalData = new ArrayList<>();

        //open file
        File file = new File(filename);
        try(BufferedReader br = new BufferedReader(new FileReader(file))) {
            String line;
            //if we need to process the cell lines
            if(isPhospho) {
                //the first line is a header
                if((line = br.readLine()) != null) {
                    String[] splitted = line.split("\t");
                    //create the cell line array
                    cellLines = new String[cellColNums.length];
                    for(int i = 0; i < cellColNums.length; i++) {
                        cellLines[i] = splitted[cellColNums[i]];
                    }
                }
            }
            //process the data
            while((line = br.readLine()) != null) {
                finalData.add(line.split("\t"));
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
     * Deletes all output files if they exist to avoid file corruption
     */
    private static void deleteAllFiles() {
        //TODO: ADD IN DELETION OF ALL OUTPUT FILES
    }

    /**
     * Deletes an individual file if it exists
     * @param filename The name of the file to delete
     */
    private static void deleteFile(String filename) {
        //check file exists
        File file = new File(filename);
        boolean exists = file.exists();

        //delete file
        if(exists) {
            file.delete();
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
            HashMap<String, HashSet<Protein>> uniProtMap
                    = new HashMap<String, HashSet<Protein>>();
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
                    HashSet<Protein> newSet = new HashSet<Protein>();
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

}

