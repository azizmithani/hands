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
import java.util.TreeSet;
import java.util.logging.Level;
import java.util.logging.Logger;

/**
 *
 * @author Aziz Mithani <aziz.mithani@lums.edu.pk>
 */
public class SNP implements Comparable<Object> {

    private int position;
    private char base;
    private int readCount;
    private static final char FIELDS_SEPERATOR = '\t';
    private static final int DIPLOID_COVERAGE_THRESHOLD = 3;
    private static final int POLYPLOID_BASE_COVERAGE_THRESHOLD = 2;
    private static final int DIPLOID_BASE_COVERAGE_THRESHOLD = 3;
    private static final double POLYPLOID_BASE_PROPORTION_THRESHOLD = 0.05;
    private static final double DIPLOID_BASE_PROPORTION_THRESHOLD = 0.30;

    public SNP() {
        position = -1;
        base = ' ';
        readCount = 0;
    }

    public SNP(int position, char base) {
        this.position = position;
        this.base = base;
        this.readCount = 0;
    }

    public SNP(int position, char base, int readCount) {
        this.position = position;
        this.base = base;
        this.readCount = readCount;
    }

    public char getBase() {
        return base;
    }

    public void setBase(char base) {
        this.base = base;
    }

    public int getPosition() {
        return position;
    }

    public void setPosition(int position) {
        this.position = position;
    }

    public int getReadCount() {
        return readCount;
    }

    public void setReadCount(int readCount) {
        this.readCount = readCount;
    }

    @Override
    public int compareTo(Object object) {
        SNP snp = (SNP) object;
        if (this.position < snp.position) {
            return -1;
        } else if (this.position > snp.position) {
            return 1;
        } else {
            if (this.base < snp.base) {
                return -1;
            } else if (this.base > snp.base) {
                return 1;
            } else {
                return 0;
            }
        }

    }

    @Override
    public int hashCode() {
        int hash = 3;
        hash = 59 * hash + this.position;
        hash = 59 * hash + this.base;
        return hash;
    }

    @Override
    public boolean equals(Object object) {
        if (!(object instanceof SNP)) {
            return false;
        }
        SNP snp = (SNP) object;
        return this.position == snp.position && this.base == snp.base;
    }

    @Override
    public String toString() {
        return position + ":" + base;// + "(" + readCount + ")";
    }

    public static HashMap initialiseBaseList(HashMap<String, Integer> sequenceLengths) {
        HashMap baseList = new HashMap();
        Iterator<Map.Entry<String, Integer>> it = sequenceLengths.entrySet().iterator();
        while (it.hasNext()) {
            Map.Entry<String, Integer> entry = it.next();

            // Get sequence details
            String seqName = entry.getKey();
            int length = entry.getValue();

            char[] seqBaseList = new char[length];
            for (int i = 0; i < length; i++) {
                seqBaseList[i] = DNA.BASE_BLANK;
            }
            baseList.put(seqName, seqBaseList);
        }

        return baseList;
    }

//    public static int readSNPList(String snpFile, HashMap sequenceLengths, HashMap snpList) {
    public static int readSNPList(String snpFile, HashMap snpList, HashMap referenceBaseList) {

        if (snpFile == null || snpFile.isEmpty()) {
            return 0;
        }

        // snpList = initialiseBaseList(sequenceLengths);
        int count = 0;
        try {
            BufferedReader br = new BufferedReader(new FileReader(snpFile));
            String previousSequenceName = "";
            char[] seqSNPList = null;
            char[] seqReferenceBaseList = null;
            String line;
            while ((line = br.readLine()) != null) {
                // ignore empty lines
                if (line.isEmpty()) {
                    continue;
                }

                int pos = 0, end;
                // sequence name 
                end = line.indexOf(FIELDS_SEPERATOR, pos);
                String sequenceName = line.substring(pos, end);
                // position
                pos = end + 1;
                end = line.indexOf(FIELDS_SEPERATOR, pos);
                int position = Integer.parseInt(line.substring(pos, end));
                // reference base
                pos = end + 1;
                end = line.indexOf(FIELDS_SEPERATOR, pos);
                String referenceBase = line.substring(pos, end);
                // consensus (SNP) base
                pos = end + 1;
                String consensusBase = line.substring(pos);

                if (!sequenceName.equals(previousSequenceName)) {
                    if (!previousSequenceName.isEmpty()) {
                        snpList.put(previousSequenceName, seqSNPList);
                        referenceBaseList.put(previousSequenceName, seqReferenceBaseList);
                    }
                    seqSNPList = (char[]) snpList.get(sequenceName);
                    seqReferenceBaseList = (char[]) referenceBaseList.get(sequenceName);
                    previousSequenceName = sequenceName;
                }

                // set the reference base
                seqReferenceBaseList[position - 1] = referenceBase.charAt(0);

                // set the base
                seqSNPList[position - 1] = consensusBase.charAt(0);
                count++;

            } // end while
            // end while
            snpList.put(previousSequenceName, seqSNPList);
        } catch (FileNotFoundException ex) {
            Logger.getLogger(SNP.class.getName()).log(Level.SEVERE, null, ex);
        } catch (IOException ex) {
            Logger.getLogger(SNP.class.getName()).log(Level.SEVERE, null, ex);
        }

        return count;
    }

    public static HashMap readReferenceBaseList(String snpFile, HashMap sequenceLengths) {
        HashMap baseList = initialiseBaseList(sequenceLengths);
        try {
            BufferedReader br = new BufferedReader(new FileReader(snpFile));
            String previousSequenceName = "";
            char[] seqBaseList = null;
            String line;
            while ((line = br.readLine()) != null) {
                // ignore empty lines
                if (line.isEmpty()) {
                    continue;
                }

                int pos = 0, end;
                // sequence name 
                end = line.indexOf(FIELDS_SEPERATOR, pos);
                String sequenceName = line.substring(pos, end);
                // position
                pos = end + 1;
                end = line.indexOf(FIELDS_SEPERATOR, pos);
                int position = Integer.parseInt(line.substring(pos, end));
                // reference base
                pos = end + 1;
                end = line.indexOf(FIELDS_SEPERATOR, pos);
                String referenceBase = line.substring(pos, end);

                if (!sequenceName.equals(previousSequenceName)) {
                    if (!previousSequenceName.isEmpty()) {
                        baseList.put(previousSequenceName, seqBaseList);
                    }
                    seqBaseList = (char[]) baseList.get(sequenceName);
                    previousSequenceName = sequenceName;
                }

                // set the base
                seqBaseList[position - 1] = referenceBase.charAt(0);
            } // end while
            // end while
            baseList.put(previousSequenceName, seqBaseList);
        } catch (FileNotFoundException ex) {
            Logger.getLogger(SNP.class.getName()).log(Level.SEVERE, null, ex);
        } catch (IOException ex) {
            Logger.getLogger(SNP.class.getName()).log(Level.SEVERE, null, ex);
        }

        return baseList;

    }

    public static int checkPolyploidHSPList(HashMap<String, char[]> polyploidHSPList,
            HashMap<String, BaseDistribution[]> polyploidBaseDistributionList, HashMap<String, char[]> referenceBaseList) {

        int count = 0;
        Iterator<Map.Entry<String, char[]>> itPolyploidHSP = polyploidHSPList.entrySet().iterator();
        while (itPolyploidHSP.hasNext()) {
            // get the next entry
            Map.Entry<String, char[]> entry = itPolyploidHSP.next();

            // get sequence name
            String sequenceName = entry.getKey();
            // get snp list for this sequence
            char[] sequenceHSPList = entry.getValue();

            // and the reference base list for this sequence
            char[] seqReferenceBaseList = referenceBaseList.get(sequenceName);
            // and the base distribution list
            BaseDistribution[] seqBaseDistribution = polyploidBaseDistributionList.get(sequenceName);

            for (int position = 0; position < sequenceHSPList.length; position++) {
                if (seqReferenceBaseList[position] == DNA.BASE_BLANK) { // no reference base assigned => no hsp or snp at this positon in any of the genome, move to the next position
                    continue;
                }
                
                // get the hsp base at this position
                char base = sequenceHSPList[position];

                // get the base distribution for this position
                BaseDistribution baseDistribution = seqBaseDistribution[position];
                int coverage = baseDistribution.getCoverage();
                
                TreeSet newBases = new TreeSet();
                if (checkPolyploidBase(baseDistribution.getCoverageA(), coverage)) {
                    newBases.add('A');
                }
                if (checkPolyploidBase(baseDistribution.getCoverageC(), coverage)) {
                    newBases.add('C');
                }
                if (checkPolyploidBase(baseDistribution.getCoverageG(), coverage)) {
                    newBases.add('G');
                }
                if (checkPolyploidBase(baseDistribution.getCoverageT(), coverage)) {
                    newBases.add('T');
                }

                if (!newBases.isEmpty()) {
                    char newBase = (Character) DNA.getReverseNucleotideMap().get(newBases.toString());
                    if (base != newBase) {
                        sequenceHSPList[position] = newBase;
                    }
                    if (DNA.isExtendedNucleotide(newBase)) {
                        count++;
                    }
                } else { // no valid base found after checking so remove the base
                    //System.out.println(position);
                    sequenceHSPList[position] = DNA.BASE_BLANK;
                }
            } // end for each position

        }
        return count;
    }

    public static void checkDiploidSNPList(HashMap<String, char[]> diploidSNPList,
            HashMap<String, BaseDistribution[]> diploidBaseDistributionList, HashMap<String, char[]> polyploidHSPList,
            HashMap<String, char[]> referenceBaseList) {

        Iterator<Map.Entry<String, char[]>> itHSPList = polyploidHSPList.entrySet().iterator();
        while (itHSPList.hasNext()) {
            // get the next entry
            Map.Entry<String, char[]> entry = itHSPList.next();

            // get sequence name
            String seqName = entry.getKey();
            // get snp list for this sequence
            char[] seqPolyploidHSPList = entry.getValue();

            // also get diploid snp list for this sequence
            char[] seqDiploidSNPList = diploidSNPList.get(seqName);
            // and the reference base list for this sequence
            char[] seqReferenceBaseList = referenceBaseList.get(seqName);
            // and the base distribution list
            BaseDistribution[] seqBaseDistribution = diploidBaseDistributionList.get(seqName);

            for (int position = 0; position < seqPolyploidHSPList.length; position++) {

                if (seqPolyploidHSPList[position] == DNA.BASE_BLANK) { // no snp in polyploid, move to the next position
                    continue;
                }

                // get the base in diploid
                char base = seqDiploidSNPList[position];
                if (base != DNA.BASE_BLANK) { // SNP has been called at this position, move to the next position
                    // check if its an ambiguous base
                    if (DNA.isExtendedNucleotide(base)) {
                        seqDiploidSNPList[position] = DNA.BASE_AMBIGUOUS;
                    }
                    //move to the next position
                    continue;
                }

                // we are now checking the position which is not called a SNP in the diploid.
                // check if multiple bases are present at this position (possibly, heterozygous)
                // get the base distribution for this position
                BaseDistribution baseDistribution = seqBaseDistribution[position];
                int coverage = baseDistribution.getCoverage();

                if (coverage == 0) {
                    seqDiploidSNPList[position] = DNA.BASE_ZERO_COVERAGE;
                    //seqDiploidSNPList[position] = seqReferenceBaseList[position];
                    continue;
                } else if (coverage < DIPLOID_COVERAGE_THRESHOLD) {
                    seqDiploidSNPList[position] = DNA.BASE_LOW_COVERAGE;
                    //seqDiploidSNPList[position] = seqReferenceBaseList[position];
                    continue;
                }

                TreeSet newBases = new TreeSet();
                if (checkDiploidBase(baseDistribution.getCoverageA(), coverage)) {
                    newBases.add('A');
                }
                if (checkDiploidBase(baseDistribution.getCoverageC(), coverage)) {
                    newBases.add('C');
                }
                if (checkDiploidBase(baseDistribution.getCoverageG(), coverage)) {
                    newBases.add('G');
                }
                if (checkDiploidBase(baseDistribution.getCoverageT(), coverage)) {
                    newBases.add('T');
                }

                if (newBases.isEmpty()) {// ideally, this should not happen 
                    seqDiploidSNPList[position] = seqReferenceBaseList[position];
                } else if (newBases.size() == 1) {
                    // add it to the list
                    seqDiploidSNPList[position] = (Character) newBases.first();
                } else if (newBases.size() > 1) { // multiple bases found at this position, mark this position as '*'
                    seqDiploidSNPList[position] = DNA.BASE_AMBIGUOUS;
                    //System.out.println(position + "\t" + base + "\t" + newBases);

                }

            } // end for each position

        }
    }

    private static boolean checkPolyploidBase(int baseCoverage, int totalCoverage) {
        if (baseCoverage >= POLYPLOID_BASE_COVERAGE_THRESHOLD && (double) baseCoverage / (double) totalCoverage >= POLYPLOID_BASE_PROPORTION_THRESHOLD) {
            return true;
        } else {
            return false;
        }

    }

    private static boolean checkDiploidBase(int baseCoverage, int totalCoverage) {
        if (baseCoverage >= DIPLOID_BASE_COVERAGE_THRESHOLD && (double) baseCoverage / (double) totalCoverage >= DIPLOID_BASE_PROPORTION_THRESHOLD) {
            return true;
        } else {
            return false;
        }

    }
}
