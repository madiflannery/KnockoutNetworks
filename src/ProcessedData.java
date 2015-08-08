import org.biopax.paxtools.model.level3.PhysicalEntity;
import org.biopax.paxtools.model.level3.Protein;

import java.util.ArrayList;
import java.util.HashMap;

/**
 * Created by Madison on 8/08/2015.
 */
public class ProcessedData {
    private ArrayList<String> geneList;
    private HashMap<String, Protein> proteins;
    private ArrayList<ArrayList<PhysicalEntity>> complexes;

    public ProcessedData() {
        this.geneList = null;
        this.proteins = null;
        this.complexes = null;
    }

    public ArrayList<String> getGeneList() {
        return geneList;
    }

    public void setGeneList(ArrayList<String> geneList) {
        this.geneList = geneList;
    }

    public HashMap<String, Protein> getProteins() {
        return proteins;
    }

    public void setProteins(HashMap<String, Protein> proteins) {
        this.proteins = proteins;
    }

    public ArrayList<ArrayList<PhysicalEntity>> getComplexes() {
        return complexes;
    }

    public void setComplexes(ArrayList<ArrayList<PhysicalEntity>> complexes) {
        this.complexes = complexes;
    }
}
