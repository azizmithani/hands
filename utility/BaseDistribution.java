/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package utility;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.HashMap;
import java.util.Iterator;
import java.util.Map;
import java.util.logging.Level;
import java.util.logging.Logger;

/**
 *
 * @author Aziz Mithani <aziz.mithani@lums.edu.pk>
 */
public class BaseDistribution {
    // constants
    private static final char FIELDS_SEPERATOR = '\t';
    // variables
    private int coverageA;
    private int coverageC;
    private int coverageG;
    private int coverageT;
    private int coverageN;

    public BaseDistribution() {
        this.coverageA = 0;
        this.coverageC = 0;
        this.coverageG = 0;
        this.coverageT = 0;
        this.coverageN = 0;
    }

    public BaseDistribution(int coverageA, int coverageT, int coverageG, int coverageC, int coverageN) {
        this.coverageA = coverageA;
        this.coverageC = coverageC;
        this.coverageG = coverageG;
        this.coverageT = coverageT;
        this.coverageN = coverageN;
    }

    public int getCoverage() {
        return getCoverageA() + getCoverageC() + getCoverageG() + getCoverageT() + getCoverageN();
    }

    public void setCoverage(int coverageA, int coverageT, int coverageG, int coverageC, int coverageN) {
        this.coverageA = coverageA;
        this.coverageC = coverageC;
        this.coverageG = coverageG;
        this.coverageT = coverageT;
        this.coverageN = coverageN;
    }
    
    public static HashMap<String, BaseDistribution[]> read(String baseDistributionFile, HashMap<String, Integer> sequenceLengths) {

        HashMap<String, BaseDistribution[]> baseDistributionList = initialise(sequenceLengths);

        try {
            BufferedReader br = new BufferedReader(new FileReader(baseDistributionFile));

            String previousSeqName = "";
            BaseDistribution[] seqBaseDistributionList = null;
            String line;
            while ((line = br.readLine()) != null) {
                // ignore empty lines
                if (line.isEmpty()) {
                    continue;
                } else if (line.startsWith("@")) { // ignore header
                    continue;
                } else {

                    int pos = 0, end;
                    // sequence name 
                    end = line.indexOf(FIELDS_SEPERATOR, pos);
                    String sequenceName = line.substring(pos, end);
                    // position
                    pos = end + 1;
                    end = line.indexOf(FIELDS_SEPERATOR, pos);
                    Integer position = Integer.parseInt(line.substring(pos, end));
                    // coverage A
                    pos = end + 1;
                    end = line.indexOf(FIELDS_SEPERATOR, pos);
                    Integer coverageA = Integer.parseInt(line.substring(pos, end));
                    // coverage T
                    pos = end + 1;
                    end = line.indexOf(FIELDS_SEPERATOR, pos);
                    Integer coverageT = Integer.parseInt(line.substring(pos, end));
                    // coverage G
                    pos = end + 1;
                    end = line.indexOf(FIELDS_SEPERATOR, pos);
                    Integer coverageG = Integer.parseInt(line.substring(pos, end));
                    // coverage C
                    pos = end + 1;
                    end = line.indexOf(FIELDS_SEPERATOR, pos);
                    Integer coverageC = Integer.parseInt(line.substring(pos, end));
                    // coverage N
                    pos = end + 1;
                    Integer coverageN = Integer.parseInt(line.substring(pos));

                    if (!sequenceName.equals(previousSeqName)) {
                        if (!previousSeqName.isEmpty()) {
                            baseDistributionList.put(previousSeqName, seqBaseDistributionList);
                        }
                        seqBaseDistributionList = (BaseDistribution[]) baseDistributionList.get(sequenceName);

                        // set the sequence name as previous sequence name
                        previousSeqName = sequenceName;
                    }
                    seqBaseDistributionList[position - 1] = new BaseDistribution(coverageA, coverageT, coverageG, coverageC, coverageN);
                }
            } // end while
            baseDistributionList.put(previousSeqName, seqBaseDistributionList);
        } catch (FileNotFoundException ex) {
            Logger.getLogger(SAM.class.getName()).log(Level.SEVERE, null, ex);
        } catch (IOException ex) {
            Logger.getLogger(SAM.class.getName()).log(Level.SEVERE, null, ex);
        }

        return baseDistributionList;

    }

    private static HashMap<String, BaseDistribution[]> initialise(HashMap<String, Integer> sequenceLengths) {

        HashMap<String, BaseDistribution[]> baseDistributionList = new HashMap();
        if (sequenceLengths == null) {
            return baseDistributionList;
        }

        Iterator<Map.Entry<String, Integer>> itSequenceLength = sequenceLengths.entrySet().iterator();
        while (itSequenceLength.hasNext()) {
            // get the next entry
            Map.Entry<String, Integer> entry = itSequenceLength.next();

            // get the seq name
            String sequenceName = entry.getKey();
            // and the length
            int length = entry.getValue();

            // initialise the base distibution list for this sequence
            BaseDistribution[] seqBaseDistributionList = new BaseDistribution[length];
            // add it to the list
            baseDistributionList.put(sequenceName, seqBaseDistributionList);
        }

        return baseDistributionList;
    }

    /**
     * @return the coverageA
     */
    public int getCoverageA() {
        return coverageA;
    }

    /**
     * @return the coverageC
     */
    public int getCoverageC() {
        return coverageC;
    }

    /**
     * @return the coverageG
     */
    public int getCoverageG() {
        return coverageG;
    }

    /**
     * @return the coverageT
     */
    public int getCoverageT() {
        return coverageT;
    }

    /**
     * @return the coverageN
     */
    public int getCoverageN() {
        return coverageN;
    }
}
