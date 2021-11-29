/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package utility;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.logging.Level;
import java.util.logging.Logger;

/**
 *
 * @author Aziz Mithani <aziz.mithani@lums.edu.pk>
 */
public class GFF {

    // constants
    public static final int GENE_FLANKING_REGION_LENGTH = 20;
    private static final String FEATURE_EXON = "exon";
    static final char FIELDS_SEPERATOR = '\t';
    // variables
    private HashMap GFFList = new HashMap();

    public GFF(String filename) {
        readGFFFile(filename);
    }

    private void readGFFFile(String filename) {

        // initialise the GFF List
        GFFList.clear();

        BufferedReader br = null;
        try {
            br = new BufferedReader(new FileReader(filename));
            String line;
            String previousSeqName = "";
            ArrayList<GFFEntry> seqGFFList = null;
            while ((line = br.readLine()) != null) {
                // ignore empty lines
                if (line.isEmpty()) {
                    continue;
                }
                // ignore the header
                if (line.startsWith("#")) {
                    continue;
                }

                int startPos = 0, endPos;
                // reference sequence name 
                endPos = line.indexOf(FIELDS_SEPERATOR, startPos);
                String sequenceName = line.substring(startPos, endPos);
                // source
                startPos = endPos + 1;
                endPos = line.indexOf(FIELDS_SEPERATOR, startPos);
                String source = line.substring(startPos, endPos);
                // feature
                startPos = endPos + 1;
                endPos = line.indexOf(FIELDS_SEPERATOR, startPos);
                String feature = line.substring(startPos, endPos);
                // start
                startPos = endPos + 1;
                endPos = line.indexOf(FIELDS_SEPERATOR, startPos);
                int start = Integer.parseInt(line.substring(startPos, endPos));
                // end
                startPos = endPos + 1;
                endPos = line.indexOf(FIELDS_SEPERATOR, startPos);
                int end = Integer.parseInt(line.substring(startPos, endPos));
                // score
                startPos = endPos + 1;
                endPos = line.indexOf(FIELDS_SEPERATOR, startPos);
                String score = line.substring(startPos, endPos);
                // strand
                startPos = endPos + 1;
                endPos = line.indexOf(FIELDS_SEPERATOR, startPos);
                char strand = line.substring(startPos, endPos).charAt(0);
                // end
                startPos = endPos + 1;
                endPos = line.indexOf(FIELDS_SEPERATOR, startPos);
                String frame = line.substring(startPos, endPos);
                // attribute
                startPos = endPos + 1;
                String attribute = line.substring(startPos);

                // save the gff entries for previous sequence if we are moving to new reference sequence
                if (!sequenceName.equals(previousSeqName)) {
                    if (!previousSeqName.isEmpty()) {
                        GFFList.put(previousSeqName, seqGFFList);
                    }

                    // initialise the GFF list for the new sequence
                    seqGFFList = new ArrayList();
                    // save the sequence name
                    previousSeqName = sequenceName;
                }

                //System.out.println(sequenceName + "\t" + Integer.toString(start) + "\t" + Integer.toString(end));
                // add it to the list
                seqGFFList.add(new GFFEntry(sequenceName, start, end, attribute, feature, strand));
                //System.out.println(line);
            }
            // add the last sequence
            GFFList.put(previousSeqName, seqGFFList);

            br.close();
        } catch (FileNotFoundException ex) {
            Logger.getLogger(GFF.class.getName()).log(Level.SEVERE, null, ex);
        } catch (IOException ex) {
            Logger.getLogger(GFF.class.getName()).log(Level.SEVERE, null, ex);
        } finally {
            try {
                br.close();
            } catch (IOException ex) {
                Logger.getLogger(GFF.class.getName()).log(Level.SEVERE, null, ex);
            }
        }
    }

    /**
     * @return the GFFList
     */
    public HashMap getGFFList() {
        return GFFList;
    }

    /**
     * @param GFFList the GFFList to set
     */
    public void setGFFList(HashMap GFFList) {
        this.GFFList = GFFList;
    }
}
