/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package hspcharacterisation;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.LinkedHashMap;
import java.util.LinkedHashSet;
import java.util.LinkedList;
import java.util.Map;
import java.util.Set;
import java.util.TreeMap;
import java.util.TreeSet;
import utility.BaseDistribution;
import utility.DNA;
import utility.GFF;
import utility.GFFEntry;
import utility.SAM;
import utility.SAMAlignment;
import utility.SNP;

/**
 *
 * @author Aziz Mithani <aziz.mithani@lums.edu.pk>
 */
public class HSPCharacterisation {

    // constant
    static final int MAX_DIPLOIDS = 10;
    // Global variables
    static double SNP_PAIR_PROPORTION_THRESHOLD = 0.05;
    static double SNP_PATTERN_MATCHING_THRESHOLD = 0.5;
    static double SINGLE_POSITION_SNP_PATTERN_COUNT_THRESHOLD = 3;
    static boolean RECTIFY_USING_REFERENCE = false;
    static int DISTANT_GENOME = 0;
    static int DIPLOID_COUNT = 0;
    static int MISSING_DIPLOID = 0;

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {

        String SAMFile = "";
        String GFFFile = "";
        String polyploidHSPFile = "";
        ArrayList<String> diploidSNPFiles = new ArrayList<String>(MAX_DIPLOIDS);
        String polyploidBaseDistributionFile = "";
        ArrayList<String> diploidBaseDistributionFiles = new ArrayList<String>(MAX_DIPLOIDS);
        ArrayList<String> outFiles = new ArrayList<String>(MAX_DIPLOIDS);

        for (int i = 0; i < MAX_DIPLOIDS; i++) {
            diploidSNPFiles.add(null);
            diploidBaseDistributionFiles.add(null);
            outFiles.add(null);
        }

        //args = "-s ../example/polyploid.sam -g ../example/reference.fa.gff -hsp ../example/polyploid.hsp -snp1 ../example/diploid1.snp -snp2 ../example/diploid2.snp -snp3 ../example/diploid3.snp -bd ../example/polyploid.bd -bd1 ../example/diploid1.bd -bd2 ../example/diploid2.bd -bd3 ../example/diploid3.bd -out1 ../example/out1.txt -out2 ../example/out2.txt -out3 ../example/out3.txt".split(" +");
        //args = "-s ../example/polyploid.sam -g ../example/reference.fa.gff -hsp ../example/polyploid.hsp -snp1 ../example/diploid1.snp -snp2 \"\" -snp3 ../example/diploid3.snp -bd ../example/polyploid.bd -bd1 ../example/diploid1.bd  -bd3 ../example/diploid3.bd -out1 ../example/out1.txt -out2 ../example/out2.txt -out3 ../example/out3.txt".split(" +");
        // read input arguments
        int i = 0;
        int nArgs = args.length;
        if (nArgs == 0) {
            printHelp();
            return;
        }
        while (i < nArgs) {
            String parameter = args[i++];
            if (parameter.equals("-s")) { // SAM file
                SAMFile = args[i++];
            } else if (parameter.equals("-g")) { // GFF file
                GFFFile = args[i++];
            } else if (parameter.equals("-hsp")) { // polyploid HSP File
                polyploidHSPFile = args[i++];
            } else if (parameter.startsWith("-snp")) { // Diploid SNP Files
                int diploidNo = Integer.parseInt(parameter.substring(4)) - 1;
                if (diploidNo < 0) {
                    System.err.println("Invalid Diploid Number.");
                    return;
                }
                String filename = args[i++];
                filename = filename.equals("\"\"") ? "" : filename;
                diploidSNPFiles.set(diploidNo, filename);
            } else if (parameter.equals("-bd")) { // polyploid Base Distribution File
                polyploidBaseDistributionFile = args[i++];
            } else if (parameter.startsWith("-bd")) { // Diploid Base Distribution Files
                int diploidNo = Integer.parseInt(parameter.substring(3)) - 1;
                if (diploidNo < 0) {
                    System.err.println("Invalid Diploid Number.");
                    return;
                }
                String filename = args[i++];
                filename = filename.equals("\"\"") ? "" : filename;
                diploidBaseDistributionFiles.set(diploidNo, filename);
            } else if (parameter.startsWith("-out")) { // Diploid Output Files
                int diploidNo = Integer.parseInt(parameter.substring(4)) - 1;
                if (diploidNo < 0) {
                    System.err.println("Invalid Diploid Number.");
                    return;
                }
                String filename = args[i++];
                filename = filename.equals("\"\"") ? "" : filename;
                outFiles.set(diploidNo, filename);
            } else if (parameter.equals("-sp")) { // SNP Pair Proportion threshold
                SNP_PAIR_PROPORTION_THRESHOLD = Double.valueOf(args[i++]);
            } else if (parameter.equals("-pm")) { // SNP Pattern Matching threshold
                SNP_PATTERN_MATCHING_THRESHOLD = Double.valueOf(args[i++]);
            } else if (parameter.equals("-r")) { // rectify using reference genome
                char rectify = args[i++].charAt(0);
                RECTIFY_USING_REFERENCE = (Character.toUpperCase(rectify) == 'T');
            } else if (parameter.equals("-d")) { //  distant genome number
                DISTANT_GENOME = Integer.valueOf(args[i++]);
            } else if (parameter.equals("-h") || parameter.equals("-help")) { // Help
                if (nArgs == 1) {
                    printHelp();
                } else {
                    System.err.println("Invalid option used with Help!");
                }
                return;
            } else {
                System.err.println("Invalid option: " + args[i - 1]);
                printHelp();
                return;
            }
        }

        // validate the input 
        if (SAMFile.isEmpty()) {
            System.err.println("SAM file must be specified.");
            return;
        } else if (GFFFile.isEmpty()) {
            System.err.println("GFF file must be specified.");
            return;
        } else if (polyploidHSPFile.isEmpty()) {
            System.err.println("Polyploid HSP File must be specified.");
            return;
        }

        // validate diploid SNP files
        for (int d = MAX_DIPLOIDS - 1; d >= 0; d--) {
            if (diploidSNPFiles.get(d) != null) {
                DIPLOID_COUNT = d + 1;
                break;
            }
        }
        if (DIPLOID_COUNT == 0) {
            System.err.println("At least one diploid SNP file must be specified.");
            return;
        }

        int missingDiploidCount = 0;
        for (int d = 0; d < DIPLOID_COUNT; d++) {
            if (diploidSNPFiles.get(d) == null) {
                System.err.println("Diploid " + Integer.toString(d + 1) + " SNP file missing. Use \"\" as the file name if you would like to skip the file.");
                return;
            } else if (diploidSNPFiles.get(d).isEmpty()) {
                missingDiploidCount++;
                MISSING_DIPLOID = d + 1;
            }
        }
        if (missingDiploidCount > 1) {
            System.err.println("Only one diploid SNP file can be skipped.");
            return;
        }
        // shrink the array list
        diploidSNPFiles.subList(DIPLOID_COUNT, diploidSNPFiles.size()).clear();

        if (DISTANT_GENOME > DIPLOID_COUNT) {
            System.err.println("Invalid Distant Genome.");
            return;
        }
        if (DISTANT_GENOME > 0 && MISSING_DIPLOID > 0) {
            System.err.println("Distant Genome cannot be used when a genome is missing.");
            return;
        } else if (MISSING_DIPLOID > 0) {
            DISTANT_GENOME = MISSING_DIPLOID;
        }

        // validate the base distribution and output files       
        for (int d = MAX_DIPLOIDS - 1; d >= 0; d--) {
            if (diploidBaseDistributionFiles.get(d) != null) {
                if (d > DIPLOID_COUNT - 1) {
                    System.err.println("Extra diploid base distribution file(s) supplied.");
                    return;
                } else if (!diploidBaseDistributionFiles.get(d).isEmpty() && d == MISSING_DIPLOID - 1) {
                    System.err.println("Diploid base distribution file for skipped SNP file supplied.");
                    return;
                }
            }
            if (outFiles.get(d) != null && d > DIPLOID_COUNT - 1) {
                System.err.println("Extra output file(s) supplied.");
                return;
            }
        }
        // shrink the array lists
        diploidBaseDistributionFiles.subList(DIPLOID_COUNT, diploidBaseDistributionFiles.size()).clear();
        outFiles.subList(DIPLOID_COUNT, outFiles.size()).clear();

        // replace null with empty string for base distribution and a suitable name for output file
        for (int d = 0; d < DIPLOID_COUNT; d++) {
            if (diploidBaseDistributionFiles.get(d) == null) {
                diploidBaseDistributionFiles.set(d, "");
            }
            if (outFiles.get(d) == null || outFiles.get(d).isEmpty()) {
                if (d + 1 == MISSING_DIPLOID) {
                    outFiles.set(d, "diploid" + Integer.toString(d + 1) + ".out");
                } else {
                    outFiles.set(d, diploidSNPFiles.get(d) + ".out");
                }
            }
        }
        int diploidBaseDistributionFilesCount = 0;
        for (int d = 0; d < DIPLOID_COUNT; d++) {
            if (!diploidBaseDistributionFiles.get(d).isEmpty()) {
                diploidBaseDistributionFilesCount++;
            }
        }
        System.out.println("Number of diploids: " + Integer.toString(DIPLOID_COUNT));
        if (MISSING_DIPLOID > 0) {
            System.out.println("Treating diploid " + Integer.toString(MISSING_DIPLOID) + " as missing diploid.");
        }

        SAM theSAM = new SAM(SAMFile);
        GFF theGFF = new GFF(GFFFile);

        // read reference base list
        //System.out.println("Reading Reference Bases ...");
        //HashMap<String, char[]> referenceBaseList = SNP.readReferenceBaseList(polyploidHSPFile, theSAM.getSequenceLengths());
        // read SNP lists
        System.out.print("Reading HSP list ... ");
        HashMap<String, char[]> referenceBaseList = SNP.initialiseBaseList(theSAM.getSequenceLengths());
        HashMap<String, char[]> polyploidHSPList = SNP.initialiseBaseList(theSAM.getSequenceLengths());
        //int polyploidHSPListSize = SNP.readSNPList(polyploidHSPFile, theSAM.getSequenceLengths(), polyploidHSPList);
        int polyploidHSPListSize = SNP.readSNPList(polyploidHSPFile, polyploidHSPList, referenceBaseList);
        System.out.println(Integer.toString(polyploidHSPListSize) + " positions.");
        if (DIPLOID_COUNT - missingDiploidCount > 1) {
            System.out.println("Reading Diploid SNP lists ...");
        } else {
            System.out.println("Reading Diploid SNP list ...");
        }
        ArrayList<HashMap<String, char[]>> diploidSNPLists = new ArrayList<HashMap<String, char[]>>(DIPLOID_COUNT);
        for (int d = 0; d < DIPLOID_COUNT; d++) {
            HashMap<String, char[]> diploidSNPList = SNP.initialiseBaseList(theSAM.getSequenceLengths());;
            //int diploidSNPListSize = SNP.readSNPList(diploidSNPFiles.get(d), theSAM.getSequenceLengths(), diploidSNPList);            
            int diploidSNPListSize = SNP.readSNPList(diploidSNPFiles.get(d), diploidSNPList, referenceBaseList);
            if (d != (MISSING_DIPLOID - 1)) {
                System.out.println("\tDiploid " + Integer.toString(d + 1) + ": " + Integer.toString(diploidSNPListSize) + " positions.");
            }
            diploidSNPLists.add(diploidSNPList);
        }

        if (polyploidBaseDistributionFile.isEmpty() && diploidBaseDistributionFilesCount == 0) {
            // skip pre-processing 
            System.out.println("Base distibution files not specified. Skipping data pre-processing ...");
        } else {
            System.out.println("Data pre-processing ...");
            // Polyploid
            if (!polyploidBaseDistributionFile.isEmpty()) {
                System.out.println("Checking HSPs ...");
                // read base distribution files
                HashMap<String, BaseDistribution[]> polyploidBaseDistributionList = BaseDistribution.read(polyploidBaseDistributionFile, theSAM.getSequenceLengths());
                // check the polyploid SNP list to see if consensus bases have been called correctly
                polyploidHSPListSize = SNP.checkPolyploidHSPList(polyploidHSPList, polyploidBaseDistributionList, referenceBaseList);
                System.out.println("\t Number of positions: " + Integer.toString(polyploidHSPListSize));
            }
            // Diploids
            for (int d = 0; d < DIPLOID_COUNT; d++) {
                String diploidBaseDistributionFile = diploidBaseDistributionFiles.get(d);
                if (!diploidBaseDistributionFile.isEmpty()) {
                    System.out.println("Checking Diploid " + Integer.toString(d + 1) + " SNPs ...");
                    // read base distribution files
                    HashMap<String, BaseDistribution[]> diploidBaseDistributionList = BaseDistribution.read(diploidBaseDistributionFile, theSAM.getSequenceLengths());
                    // check the diploid SNP list to see if any SNP is missing 
                    SNP.checkDiploidSNPList(diploidSNPLists.get(d), diploidBaseDistributionList, polyploidHSPList, referenceBaseList);
                }
            }
        }
        ArrayList<BufferedWriter> outList = new ArrayList<BufferedWriter>(DIPLOID_COUNT);
        try {
            for (int d = 0; d < DIPLOID_COUNT; d++) {
                //Construct the BufferedWriter object
                outList.add(new BufferedWriter(new FileWriter(outFiles.get(d))));
            }
        } catch (IOException ex) {
            System.err.println("Unable to open output file(s)");
            return;
        }

        System.out.println("Processing SAM file ...");
        // process SAM file
        processSAMFile(theSAM, theGFF, polyploidHSPList, diploidSNPLists, referenceBaseList, outList);

        try {
            for (int d = 0; d < DIPLOID_COUNT; d++) {
                //close the BufferedWriter object
                outList.get(d).close();
            }
        } catch (IOException ex) {
            System.err.println("Unable to close output file(s)");
        }
    }

    private static void processSAMFile(SAM theSAM, GFF theGFF, HashMap<String, char[]> polyploidHSPList,
            ArrayList<HashMap<String, char[]>> diploidSNPLists, HashMap<String, char[]> referenceBaseList,
            ArrayList<BufferedWriter> outList) {

        SAMAlignment theSAMAlignment;
        GFFEntry theGFFEntry = null;
        ArrayList theSequenceGFFList;
        Iterator<GFFEntry> itGFFEntry = null;
        char[] seqPolyploidHSPList = null;
        char[] seqReferenceBaseList = null;
        boolean withinValidGene = false;
        LinkedHashMap<Integer, SNP> geneHSPList = null;
        LinkedHashMap<Integer, SNP> geneReferenceBaseList = null;
        HashMap<String, SAMAlignment> pendingAlignments = new HashMap();
        LinkedList<SAMAlignment> alignmentsToBeProcessed = new LinkedList();
        ArrayList<char[]> seqDiploidSNPLists = new ArrayList<char[]>(DIPLOID_COUNT);
        for (int d = 0; d < DIPLOID_COUNT; d++) {
            seqDiploidSNPLists.add(null);
        }
        ArrayList<LinkedHashMap<Integer, SNP>> geneDiploidSNPLists = new ArrayList<LinkedHashMap<Integer, SNP>>(DIPLOID_COUNT);
        for (int d = 0; d < DIPLOID_COUNT; d++) {
            geneDiploidSNPLists.add(null);
        }

        while ((theSAMAlignment = theSAM.readNextAlignment()) != null) {

            // get the GFF list for the reference against which this read is aligned
            if (theGFFEntry == null
                    || !theSAMAlignment.getReferenceSequenceName().equals(theGFFEntry.getSequenceName())
                    || theSAMAlignment.getAlignmentStart() > theGFFEntry.getEnd()
                    || theSAMAlignment.getAlignmentEnd() < theGFFEntry.getStart()) {

                if (theGFFEntry != null && withinValidGene) {
                    // process the previous set of alignments
                    processAlignments(alignmentsToBeProcessed, geneHSPList, geneDiploidSNPLists, geneReferenceBaseList, theGFFEntry.getSequenceName(), outList);
                    // reset the variables
                    alignmentsToBeProcessed.clear();
                    pendingAlignments.clear();
                    withinValidGene = false;
                }

                if (theGFFEntry == null || !theSAMAlignment.getReferenceSequenceName().equals(theGFFEntry.getSequenceName())) {

                    // get the HSP list for this sequence
                    seqPolyploidHSPList = (char[]) polyploidHSPList.get(theSAMAlignment.getReferenceSequenceName());
                    // and the SNP lists for diploids
                    for (int d = 0; d < DIPLOID_COUNT; d++) {
                        HashMap<String, char[]> diploidSNPList = diploidSNPLists.get(d);
                        if (diploidSNPList != null) {
                            seqDiploidSNPLists.set(d, (char[]) diploidSNPList.get(theSAMAlignment.getReferenceSequenceName()));
                        } else {
                            seqDiploidSNPLists.set(d, null);
                        }
                    }
                    // and the list of reference base
                    seqReferenceBaseList = (char[]) referenceBaseList.get(theSAMAlignment.getReferenceSequenceName());

                    // get the GFF list for this sequence
                    theSequenceGFFList = (ArrayList) theGFF.getGFFList().get(theSAMAlignment.getReferenceSequenceName());
                    // and the iterator to the GFF entries
                    itGFFEntry = theSequenceGFFList.iterator();
                    if (itGFFEntry.hasNext()) {
                        theGFFEntry = itGFFEntry.next();
                    } else {
                        theGFFEntry = null;
                        continue;
                    }
                }

                // move to the gene in which this read lies
                while (theSAMAlignment.getAlignmentStart() > theGFFEntry.getEnd()) {

                    if (itGFFEntry.hasNext()) {
                        theGFFEntry = itGFFEntry.next();
                    } else {
                        theGFFEntry = null;
                        break;
                    }
                }

                if (theGFFEntry == null) {
                    // this read pair lies after the last gene of this sequence. Ignore it.
                    continue;
                } else if (theSAMAlignment.getAlignmentEnd() < theGFFEntry.getStart()) {
                    // this read does not lie within a valid gene. Ignore it.
                    withinValidGene = false;
                    continue;
                }

                //System.out.println(theGFFEntry.toString());
                System.out.print(theGFFEntry.toString());

                // Now theGFFEntry corresponds to the gene in which the read lies.
                // Get the SNPs/HSPs within this gene
                geneHSPList = getGeneSNPs(seqPolyploidHSPList, theGFFEntry.getStart() - GFF.GENE_FLANKING_REGION_LENGTH, theGFFEntry.getEnd() + GFF.GENE_FLANKING_REGION_LENGTH, true);

                // also output the number of HSPs in this region
                System.out.println("\t" + Integer.toString(geneHSPList.size()));

                for (int d = 0; d < DIPLOID_COUNT; d++) {
                    char[] seqDiploidSNPList = seqDiploidSNPLists.get(d);
                    if (seqDiploidSNPList != null) {
                        geneDiploidSNPLists.set(d, getGeneSNPs(seqDiploidSNPList, theGFFEntry.getStart() - GFF.GENE_FLANKING_REGION_LENGTH, theGFFEntry.getEnd() + GFF.GENE_FLANKING_REGION_LENGTH, false));
                    } else {
                        geneDiploidSNPLists.set(d, null);
                    }
                }
                // and the reference base list for this gene
                geneReferenceBaseList = getGeneSNPs(seqReferenceBaseList, theGFFEntry.getStart() - GFF.GENE_FLANKING_REGION_LENGTH, theGFFEntry.getEnd() + GFF.GENE_FLANKING_REGION_LENGTH, false);
            }
            // reset the flag
            withinValidGene = true;

            // check if the read's mate is mapped or not
            SAMAlignment theMate;
            if ((theSAMAlignment.getFlag() & SAM.FLAG_READ_PAIRED) == 0) { // unpaired data. add to the list of alignments to be processed
                alignmentsToBeProcessed.add(theSAMAlignment);
            } else {
                if ((theSAMAlignment.getFlag() & SAM.FLAG_MATE_UNMAPPED) > 0) {
                    alignmentsToBeProcessed.add(theSAMAlignment);
                } else if ((theMate = pendingAlignments.get(theSAMAlignment.getReadName())) != null) { // see if the read's mate has been seen before
                    // Yes. Set the mate information
                    theSAMAlignment.setMate(theMate);
                    theMate.setMate(theSAMAlignment);
                    // add the mate to the list of alignments to be processed because the mate is actually mapped before this alignment
                    alignmentsToBeProcessed.add(theMate);
                    // remove the alignment from the list of pending alignments
                    pendingAlignments.remove(theSAMAlignment.getReadName());
                } else { //no, add it to the list of pending alignments
                    pendingAlignments.put(theSAMAlignment.getReadName(), theSAMAlignment);
                }
            }
        } // end while

        // process the final set of alignments 
        if (theGFFEntry != null && withinValidGene) {
            processAlignments(alignmentsToBeProcessed, geneHSPList, geneDiploidSNPLists, geneReferenceBaseList, theGFFEntry.getSequenceName(), outList);
            // reset the variables
            alignmentsToBeProcessed.clear();
            pendingAlignments.clear();
        }

    }

    private static void processAlignments(LinkedList<SAMAlignment> alignmentsToBeProcessed,
            LinkedHashMap<Integer, SNP> geneHSPList, ArrayList<LinkedHashMap<Integer, SNP>> geneDiploidSNPLists,
            LinkedHashMap<Integer, SNP> geneReferenceBaseList, String sequenceName,
            ArrayList<BufferedWriter> outList) {

        LinkedHashMap<TreeSet, Integer> basePatterns = new LinkedHashMap();
        Iterator<SAMAlignment> itSAMAlignment = alignmentsToBeProcessed.iterator();
        while (itSAMAlignment.hasNext()) {
            SAMAlignment theSAMAlignment = itSAMAlignment.next();
            SAMAlignment theMate = theSAMAlignment.getMate();

            // get the HSPs that lie within the region of these alignments
            LinkedHashSet<SNP> alignmentHSPList = getSNPs(geneHSPList, theSAMAlignment.getAlignmentStart(), theSAMAlignment.getAlignmentEnd());
            LinkedHashSet<SNP> mateHSPList = (theMate == null ? new LinkedHashSet() : getSNPs(geneHSPList, theMate.getAlignmentStart(), theMate.getAlignmentEnd()));

            // ignore the alignments if none of them have any HSPs within their boundaries
            if (alignmentHSPList.isEmpty() && mateHSPList.isEmpty()) {
                continue;
            }

            // get SNPs in the alignment read
            TreeSet alignmentBasePattern = SAM.getMismatchingBases(alignmentHSPList, theSAMAlignment.getAlignmentStart(), theSAMAlignment.getMismatchingPositions(), theSAMAlignment.getReadSequence(), theSAMAlignment.getBaseQualities(), theSAMAlignment.getCigar());
            if (alignmentBasePattern == null) { // poor data present at SNP positions .. ignore
                continue;
            }

            if (theMate != null) {
                // get SNPs in the mate read
                TreeSet mateBasePattern = SAM.getMismatchingBases(mateHSPList, theMate.getAlignmentStart(), theMate.getMismatchingPositions(), theMate.getReadSequence(), theMate.getBaseQualities(), theMate.getCigar());
                if (mateBasePattern == null) { // poor data present at SNP positions .. ignore
                    continue;
                }

                // check if there is any descrepancy in the SNPs in the two reads. This will only happen if the reads overlap
                if (theSAMAlignment.getAlignmentEnd() >= theMate.getAlignmentStart() && checkDiscrepancyInReadPair(alignmentBasePattern, mateBasePattern)) {
                    continue;
                }

                // add Read 2 SNPs to the Read 1 SNP list to get a combined SNP list
                alignmentBasePattern.addAll(mateBasePattern);
            }

            if (alignmentBasePattern.size() > 0) {
                // we have the list of valid SNPs present in this read pair.

                // Check if this snp pattern has been seen before
                Integer frequency = basePatterns.get(alignmentBasePattern);
                if (frequency == null) { // if not, add it to the list of snp patterns
                    basePatterns.put(alignmentBasePattern, 1);
                } else { // otherwise increment the number of reads containing this SNP pattern
                    basePatterns.put(alignmentBasePattern, frequency + 1);
                }
            } // if Read SNP list size > 0 

        } // end while

        // process groups for previous sequence / previous gene
        processBasePatterns(basePatterns, geneHSPList, geneDiploidSNPLists, geneReferenceBaseList, sequenceName, outList);
    }

    private static LinkedHashMap<Integer, SNP> getGeneSNPs(char[] seqSNPList, int start, int end, boolean onlyHSPs) {

        LinkedHashMap<Integer, SNP> snpList = new LinkedHashMap();

        start = Math.max(start, 1);
        end = Math.min(end, seqSNPList.length);
        for (int i = start - 1; i < end; i++) { // we start from start - 1 because SNP list is 0-based
            if (seqSNPList[i] != DNA.BASE_BLANK) {
                if (onlyHSPs && !DNA.isExtendedNucleotide(seqSNPList[i])) { // if onlyHSPs are to be returned and this is an unambiguous position then ignore it
                    continue;
                }
                snpList.put(i + 1, new SNP(i + 1, seqSNPList[i])); //add to the actual position (1-based)
            }
        }

        return snpList;
    }

    private static LinkedHashSet<SNP> getSNPs(LinkedHashMap<Integer, SNP> theSNPList, int start, int end) {

        LinkedHashSet<SNP> theSNPSubList = new LinkedHashSet();
        Iterator<Map.Entry<Integer, SNP>> itSNP = theSNPList.entrySet().iterator();
        while (itSNP.hasNext()) {
            Map.Entry<Integer, SNP> entry = itSNP.next();
            int position = entry.getKey();

            if (position < start) {
                // do nothing
            } else if (position > end) {
                break;
            } else {
                theSNPSubList.add(entry.getValue());
            }
        }

        return theSNPSubList;
    }

    private static boolean checkDiscrepancyInReadPair(TreeSet alignment1BasePattern, TreeSet alignment2BasePattern) {

        Iterator<SNP> itRead1SNPList = alignment1BasePattern.iterator();
        Iterator<SNP> itRead2SNPList = alignment2BasePattern.iterator();
        SNP snp1, snp2;
        if (itRead1SNPList.hasNext()) {
            snp1 = itRead1SNPList.next();
        } else {
            return false;
        }
        if (itRead2SNPList.hasNext()) {
            snp2 = itRead2SNPList.next();
        } else {
            return false;
        }
        while (true) {

            if (snp1.getPosition() < snp2.getPosition()) {
                // get next SNP
                if (itRead1SNPList.hasNext()) {
                    snp1 = itRead1SNPList.next();
                } else {
                    break;
                }
            } else if (snp1.getPosition() == snp2.getPosition()) { // same position
                if (snp1.getBase() == snp2.getBase()) { //same SNP
                    if (itRead1SNPList.hasNext() && itRead2SNPList.hasNext()) {
                        snp1 = itRead1SNPList.next();
                        snp2 = itRead2SNPList.next();
                    } else {
                        break;
                    }
                } else {
                    return true;
                }
            } else if (snp1.getPosition() > snp2.getPosition()) {
                if (itRead2SNPList.hasNext()) {
                    snp2 = itRead2SNPList.next();
                } else {
                    break;
                }
            }
        } // end while

        return false;
    }

    private static void processBasePatterns(LinkedHashMap<TreeSet, Integer> basePatterns,
            LinkedHashMap<Integer, SNP> geneHSPList, ArrayList<LinkedHashMap<Integer, SNP>> geneDiploidSNPLists,
            LinkedHashMap<Integer, SNP> geneReferenceBaseList, String sequenceName,
            ArrayList<BufferedWriter> outList) {

        // filter the base patterns to remove unwanted patterns
        filterBasePatterns(basePatterns, geneHSPList);

        // sort the base patterns by size (longest patterns first)
        LinkedHashSet<TreeSet> sortedBasePatterns = sortBasePatternsBySize(basePatterns);

        // remove base patterns embedded in another pattern
        removeEmbeddedBasePatterns(sortedBasePatterns);

        // assign patterns to 1st, 2nd and/or 3rd genomes
        ArrayList<LinkedHashMap<TreeSet, Double>> basePatternsList = new ArrayList<LinkedHashMap<TreeSet, Double>>(DIPLOID_COUNT);
        for (int d = 0; d < DIPLOID_COUNT; d++) {
            basePatternsList.add(new LinkedHashMap());
        }
        assignBasePatternsToGenomes(sortedBasePatterns, geneDiploidSNPLists, basePatternsList, false);

        // assign genome bases
        ArrayList<LinkedHashMap<Integer, SNP>> baseLists = new ArrayList<LinkedHashMap<Integer, SNP>>(DIPLOID_COUNT);
        for (LinkedHashMap<TreeSet, Double> genomeBasePatterns : basePatternsList) {
            baseLists.add(assignBasesToGenome(genomeBasePatterns));
        }

        // get the list of base patterns not assigned to any genome
        LinkedHashSet<TreeSet> unassignedBasePatterns = new LinkedHashSet(sortedBasePatterns);
        for (LinkedHashMap<TreeSet, Double> genomeBasePatterns : basePatternsList) {
            unassignedBasePatterns.removeAll(genomeBasePatterns.keySet());
        }

        ArrayList<Boolean> genomeFinalised = new ArrayList<Boolean>(DIPLOID_COUNT);
        for (int d = 0; d < DIPLOID_COUNT; d++) {
            genomeFinalised.add(false);
        }

        if (RECTIFY_USING_REFERENCE) { // correction using reference genome
            // check which genome(s) contributed to the reference sequence
            // we use a set as multiple genome might have the maximum identity with the reference
            HashSet referenceGenome = checkReferenceGenome(geneHSPList, geneReferenceBaseList, baseLists);

            for (int d = 0; d < DIPLOID_COUNT; d++) {
                if (referenceGenome.contains(d)) { // reference comes from this genome 
                    // rectify the SNP list for this genome using reference base list
                    rectifyBaseListUsingReferenceGenome(geneHSPList, geneReferenceBaseList, geneDiploidSNPLists.get(d), baseLists.get(d));
                    // remove the patterns which no longer contradict the base list (i.e. no longer unassigned)
                    checkBasePatternsForConsistency(unassignedBasePatterns, baseLists.get(d), true);
                    genomeFinalised.set(d, true);
                }
            }
        }

        // assign unassigned SNP patterns to genomes
        ArrayList<LinkedHashMap<TreeSet, Double>> unassignedBasePatternsList = new ArrayList<LinkedHashMap<TreeSet, Double>>(DIPLOID_COUNT);
        for (int d = 0; d < DIPLOID_COUNT; d++) {
            unassignedBasePatternsList.add(new LinkedHashMap());
        }
        assignBasePatternsToGenomes(unassignedBasePatterns, baseLists, unassignedBasePatternsList, true);

        // finalise the base lists, if need and then write the base list
        Iterator<LinkedHashMap<Integer, SNP>> itBaseLists = baseLists.iterator();
        Iterator<LinkedHashMap<TreeSet, Double>> itBasePatternsList = basePatternsList.iterator();
        Iterator<LinkedHashMap<TreeSet, Double>> itUnassignedBasePatternsList = unassignedBasePatternsList.iterator();
        Iterator<Boolean> itGenomeFinalised = genomeFinalised.iterator();
        Iterator<LinkedHashMap<Integer, SNP>> itGeneDiploidSNPLists = geneDiploidSNPLists.iterator();
        Iterator<BufferedWriter> itOutList = outList.iterator();
        while (itBaseLists.hasNext()) {
            LinkedHashMap<Integer, SNP> baseList = itBaseLists.next();
            if (!itGenomeFinalised.next()) {
                finaliseBaseList(itUnassignedBasePatternsList.next(), baseList, itBasePatternsList.next());
            }
            // commented on 2-jan-2015 1330 - writing is now done once all the processing is complete. This will allow futher processing/validation, if required.
            // write the output
            //writeBases(baseList, geneHSPList, geneReferenceBaseList, itGeneDiploidSNPLists.next(), sequenceName, itOutList.next());
        }

        LinkedHashMap<Integer, Character> hspBasesCovered = checkIfAllHSPBasesAreCovered(geneHSPList, baseLists);
        // added on 2-jan-2015 1330 - writing is now done once all the processing is complete. This will allow futher processing/validation, if required.
        // write the output for each base list
        itBaseLists = baseLists.iterator();
        while (itBaseLists.hasNext()) {
            // write the output for the next base list
            writeBases(itBaseLists.next(), geneHSPList, geneReferenceBaseList, itGeneDiploidSNPLists.next(), sequenceName, itOutList.next(), hspBasesCovered);
        }

    }

    private static void filterBasePatterns(LinkedHashMap<TreeSet, Integer> basePatterns, LinkedHashMap<Integer, SNP> geneHSPList) {

        // get SNP pairs along with the number of reads supporting each pair
        HashMap<SNPPair, Integer> snpPairList = getSNPPairs(basePatterns);
        // get SNP position pairs along with the number of reads supporting each position pair
        HashMap<PositionPair, Integer> positionPairList = getPositionPairs(basePatterns);

        HashSet basePatternsToBeDeleted = new HashSet();
        Iterator<Map.Entry<TreeSet, Integer>> itBasePatterns = basePatterns.entrySet().iterator();
        while (itBasePatterns.hasNext()) {
            // get the next entry
            Map.Entry<TreeSet, Integer> entry = itBasePatterns.next();

            // Get the next SNP pattern
            TreeSet basePattern = entry.getKey();

            // check if is a single position pattern
            if (basePattern.size() == 1) {
                // Yes! check the threshold
                if (entry.getValue() < SINGLE_POSITION_SNP_PATTERN_COUNT_THRESHOLD) {
                    basePatternsToBeDeleted.add(basePattern);
                }
                continue;
            }
            // now check the multi position base patterns
            // get the iterator over this Base pattern
            Iterator<SNP> itSNP = basePattern.iterator();
            // Get the first SNP in the Base Pattern
            SNP previousSNP = itSNP.next();
            //SNP previousSNP = null;
            //if (itSNP.hasNext()) {
            // previousSNP = itSNP.next();
            // }
            // process all remaining SNPs
            while (itSNP.hasNext()) {
                SNP theSNP = itSNP.next();

                // create a new SNP position pair using the previous and current SNPs
                PositionPair positionPair = new PositionPair(previousSNP.getPosition(), theSNP.getPosition());

                if (consecutivePositions(positionPair, geneHSPList)) {

                    // create a new SNP pair using the previous and current SNPs
                    SNPPair snpPair = new SNPPair(previousSNP, theSNP);

                    // get the readcount for this SNP pair
                    int pairReadCount = snpPairList.get(snpPair);

                    // create a new SNP position pair using the previous and current SNPs and get the read count for this snp position pair
                    int positionReadCount = positionPairList.get(positionPair);

                    if ((double) pairReadCount / (double) positionReadCount < SNP_PAIR_PROPORTION_THRESHOLD) {
                        basePatternsToBeDeleted.add(basePattern);
                        break;
                    }
                }

                // save this SNP as previous SNP
                previousSNP = theSNP;
            } // end while itSNPPattern.hasNext()

        } // it.hasNext()

        // remove the groups which have been marked for deletion
        basePatterns.keySet().removeAll(basePatternsToBeDeleted);

    }

    private static HashMap getSNPPairs(LinkedHashMap<TreeSet, Integer> basePatterns) {

        HashMap<SNPPair, Integer> snpPairList = new HashMap();
        Iterator<Map.Entry<TreeSet, Integer>> itBasePattern = basePatterns.entrySet().iterator();
        while (itBasePattern.hasNext()) {
            // Get the current SNP pattern and corresponding read counts
            Map.Entry<TreeSet, Integer> entry = itBasePattern.next();

            TreeSet basePattern = entry.getKey();
            int readCount = entry.getValue();

            SNP previousSNP = null;
            // get the iterator over this SNP pattern
            Iterator<SNP> itSNP = basePattern.iterator();
            // Get the first SNP in the SNP Pattern
            if (itSNP.hasNext()) {
                previousSNP = itSNP.next();
            }
            // process all remaining SNPs
            while (itSNP.hasNext()) {
                SNP theSNP = itSNP.next();

                // create a new SNP pair using the previous and current SNPs
                SNPPair snpPair = new SNPPair(previousSNP, theSNP);

                // get the entry from the list of snp pairs
                Integer currentReadCount = snpPairList.get(snpPair);
                if (currentReadCount == null) {
                    // add it to the list of SNP pairs
                    snpPairList.put(snpPair, readCount);
                } else {
                    // add it to the list of SNP pairs
                    snpPairList.put(snpPair, currentReadCount + readCount);
                }
                // save this SNP as previous SNP
                previousSNP = theSNP;
            }
        } // end while (snpPatterns)

        return snpPairList;
    }

    private static HashMap getPositionPairs(LinkedHashMap<TreeSet, Integer> basePatterns) {

        HashMap<PositionPair, Integer> positionPairList = new HashMap();
        Iterator<Map.Entry<TreeSet, Integer>> itBasePattern = basePatterns.entrySet().iterator();
        while (itBasePattern.hasNext()) {
            // Get the current SNP pattern and corresponding read counts
            Map.Entry<TreeSet, Integer> entry = itBasePattern.next();

            TreeSet basePattern = entry.getKey();
            int readCount = entry.getValue();

            // get the iterator over this SNP pattern
            Iterator<SNP> itSNP = basePattern.iterator();
            // Get the first SNP in the SNP Pattern
//            SNP previousSNP = null;
//            if (itSNP.hasNext()) {
//                previousSNP = itSNP.next();
//            }
            SNP previousSNP = itSNP.next();
            // process all remaining SNPs
            while (itSNP.hasNext()) {
                SNP theSNP = itSNP.next();

                // create a new SNP pair using the previous and current SNPs
                PositionPair positionPair = new PositionPair(previousSNP.getPosition(), theSNP.getPosition());

                // get the entry from the list of snp pairs
                Integer currentReadCount = positionPairList.get(positionPair);
                if (currentReadCount == null) {
                    // add it to the list of SNP pairs
                    positionPairList.put(positionPair, readCount);
                } else {
                    // add it to the list of SNP pairs
                    positionPairList.put(positionPair, currentReadCount + readCount);
                }

                // save this SNP as previous SNP
                previousSNP = theSNP;
            }
            // System.out.println();
        } // end while (snpPatterns)

        return positionPairList;
    }

    private static boolean consecutivePositions(PositionPair positionPair, LinkedHashMap<Integer, SNP> geneHSPList) {
        Iterator<Integer> itPosition = geneHSPList.keySet().iterator();
        while (itPosition.hasNext()) {
            int position = itPosition.next();
            if (positionPair.getPosition1() == position) {
                if (itPosition.hasNext()) {
                    position = itPosition.next();
                    if (positionPair.getPosition2() == position) {
                        return true;
                    } else {
                        return false;
                    }
                } else {
                    return false;
                }
            } else if (positionPair.getPosition1() < position) {
                // ideally this should not happen
                return false;
            }
        }
        return false;
    }

    private static LinkedHashSet<TreeSet> sortBasePatternsBySize(LinkedHashMap<TreeSet, Integer> basePatterns) {
        /* Sort the base patterns such that largest pattern is at the begining */

        TreeMap<Integer, LinkedHashSet> basePatternsBySize = new TreeMap();
        Iterator<TreeSet> itBasePatterns = basePatterns.keySet().iterator();
        while (itBasePatterns.hasNext()) {
            // get the next base pattern
            TreeSet basePattern = itBasePatterns.next();

            // get the list of patterns of this size
            LinkedHashSet<TreeSet> sizeBasePatterns = basePatternsBySize.get(basePattern.size());
            if (sizeBasePatterns == null) { // if no pattern found, create a new list
                sizeBasePatterns = new LinkedHashSet();
            }

            // add this pattern to the list of patterns
            sizeBasePatterns.add(basePattern);
            basePatternsBySize.put(basePattern.size(), sizeBasePatterns);
        } // end while()

        // now order them descendingly according to their size
        LinkedHashSet<TreeSet> sortedBasePatterns = new LinkedHashSet();
        Iterator<Integer> itSize = basePatternsBySize.descendingKeySet().iterator();
        while (itSize.hasNext()) {
            sortedBasePatterns.addAll(basePatternsBySize.get(itSize.next()));
        }

        return sortedBasePatterns;
    }

    private static void removeEmbeddedBasePatterns(LinkedHashSet<TreeSet> basePatterns) {

        HashSet basePatternsToBeDeleted = new HashSet();
        Iterator<TreeSet> itBasePatterns1 = basePatterns.iterator();
        while (itBasePatterns1.hasNext()) {
            // get the next SNP group
            TreeSet<SNP> basePattern1 = itBasePatterns1.next();

            // move to the next pattern if this group has been marked as inactive
            if (basePatternsToBeDeleted.contains(basePattern1)) {
                continue;
            }

            // get the first and last SNP position in this group
            int firstSNPPositionGroup1 = basePattern1.first().getPosition();
            int lastSNPPositionGroup1 = basePattern1.last().getPosition();

            Iterator<TreeSet> itBasePatterns2 = basePatterns.iterator();
            while (itBasePatterns2.hasNext()) {
                TreeSet basePattern2 = itBasePatterns2.next();
                if (basePattern1 == basePattern2) {
                    break;
                }
            }
            while (itBasePatterns2.hasNext()) {
                // get the SNP pattern
                TreeSet<SNP> basePattern2 = itBasePatterns2.next();

                if (lastSNPPositionGroup1 < basePattern2.first().getPosition()) {// if the group 1 ends before group 2 starts, no overlap possible
                    continue;
                } else if (basePattern2.last().getPosition() < firstSNPPositionGroup1) {// if the group 2 ends before group 1 starts, no overlap possible
                    continue;
                }

                int overlappingSNPs = 0;
                Iterator<SNP> itSNPPattern1 = basePattern1.iterator();
                Iterator<SNP> itSNPPattern2 = basePattern2.iterator();
//                SNP snp1 = null, snp2 = null;
//                if (itSNPPattern1.hasNext()) {
//                    snp1 = itSNPPattern1.next();
//                }
//                if (itSNPPattern2.hasNext()) {
//                    snp2 = itSNPPattern2.next();
//                }
                SNP snp1 = itSNPPattern1.next();
                SNP snp2 = itSNPPattern2.next();
                while (true) {
                    if (snp1.getPosition() == snp2.getPosition()) { // same position
                        if (snp1.getBase() == snp2.getBase()) { //same SNP
                            overlappingSNPs++;
                            // get next SNPs
                            if (itSNPPattern1.hasNext() && itSNPPattern2.hasNext()) {
                                snp1 = itSNPPattern1.next();
                                snp2 = itSNPPattern2.next();
                            } else {
                                break;
                            }
                        } else {
                            break;
                        }
                    } else if (snp1.getPosition() < snp2.getPosition()) {
                        if (itSNPPattern1.hasNext()) {
                            snp1 = itSNPPattern1.next();
                        } else {
                            break;
                        }
                    } else if (snp1.getPosition() > snp2.getPosition()) {
                        if (itSNPPattern2.hasNext()) {
                            snp2 = itSNPPattern2.next();
                        } else {
                            break;
                        }
                    }
                } // end while

                if (overlappingSNPs == basePattern2.size()) { // all SNPs are embedded in group 1
                    // add it to the groups to be deleted
                    basePatternsToBeDeleted.add(basePattern2);
                }
            } // end while itSNPGroups2.hasNext()
        }// end while itSNPGroups1.hasNext()

        // remove the groups which were found to be embedded in other groups
        basePatterns.removeAll(basePatternsToBeDeleted);
    }

    private static void assignBasePatternsToGenomes(LinkedHashSet basePatterns,
            ArrayList<LinkedHashMap<Integer, SNP>> geneDiploidSNPLists,
            ArrayList<LinkedHashMap<TreeSet, Double>> basePatternsList,
            boolean ignoreUnassingedOrAmbiguousPositions) {

        // iterate through all SNP Patterns
        Iterator<TreeSet> it = basePatterns.iterator();
        while (it.hasNext()) {
            int genomeAssignedCount = 0;
            // Get the next SNP pattern
            TreeSet basePattern = it.next();

            ArrayList<Double> proportionList = new ArrayList<Double>(DIPLOID_COUNT);
            ArrayList<Boolean> isGenomeList = new ArrayList<Boolean>(DIPLOID_COUNT);

            // check the pattern for each genome
            for (int d = 0; d < DIPLOID_COUNT; d++) {
                LinkedHashMap<Integer, SNP> geneDiploidSNPList = geneDiploidSNPLists.get(d);
                if (geneDiploidSNPList != null) {
                    double proportion = calculateBasePatternIdentity(basePattern, geneDiploidSNPList, ignoreUnassingedOrAmbiguousPositions);
                    proportionList.add(proportion);
                    // check if the matching proportion meets the threshold
                    boolean isGenome = (proportion >= SNP_PATTERN_MATCHING_THRESHOLD);
                    isGenomeList.add(isGenome);
                    if (isGenome && d != (DISTANT_GENOME - 1)) {
                        genomeAssignedCount++;
                    }
                } else {
                    proportionList.add(0.0);
                    isGenomeList.add(false);
                }
            }

            if (DISTANT_GENOME > 0 && genomeAssignedCount == 0 && basePattern.size() > 1) {
                basePatternsList.get(DISTANT_GENOME - 1).put(basePattern, 2.0);
            } else {
                for (int d = 0; d < DIPLOID_COUNT; d++) {
                    if (isGenomeList.get(d)) {
                        basePatternsList.get(d).put(basePattern, proportionList.get(d));
                    }
                }
            }
        } // end while (snpPatterns)

    }

    private static double calculateBasePatternIdentity(TreeSet basePattern, LinkedHashMap<Integer, SNP> snpList, boolean ignoreUnassingedOrAmbigousPositions) {
        int matchingSNPs = 0;
        int positionsToBeIgnored = 0;
        Iterator<SNP> itSNPPattern = basePattern.iterator();
        while (itSNPPattern.hasNext()) {
            // get the next SNP in the pattern
            SNP thePatternSNP = itSNPPattern.next();
            // get the corresponding SNP in the SNP list
            SNP theListSNP = snpList.get(thePatternSNP.getPosition());

            if (theListSNP == null || !DNA.isNucleotide(theListSNP.getBase())) {
                positionsToBeIgnored++;
                continue;
            }

            // see if the snp base in the pattern matches the snp in the list
            if (thePatternSNP.getBase() == theListSNP.getBase()) {
                // if yes, increment the match count
                matchingSNPs++;
            }
        } // end while it.hasNext()

        if (ignoreUnassingedOrAmbigousPositions) {
            return (double) matchingSNPs / (double) basePattern.size();
        } else {
            // return the match proportion
            return (double) matchingSNPs / (double) (basePattern.size() - positionsToBeIgnored);
        }

    }

    private static LinkedHashMap<Integer, SNP> assignBasesToGenome(LinkedHashMap<TreeSet, Double> potentialBasePatterns) {

        LinkedHashMap<Integer, SNP> baseList = new LinkedHashMap();

        if (potentialBasePatterns.isEmpty()) {
            return baseList;
        }

        boolean baseListChanged;
        do {

            TreeMap<Integer, HashMap> potentialBases = new TreeMap();
            // iterate through all SNP Patterns
            Iterator<Map.Entry<TreeSet, Double>> itSNPPattern = potentialBasePatterns.entrySet().iterator();
            while (itSNPPattern.hasNext()) {
                // get the next entry
                Map.Entry<TreeSet, Double> entry = itSNPPattern.next();
                // Get the next SNP pattern
                TreeSet snpPattern = entry.getKey();
                // update potential bases
                updatePotentialBases(snpPattern, entry.getValue(), potentialBases, baseList);

            } // end while (snpPatterns)

            // assign bases
            baseListChanged = assignBases(potentialBases, baseList);

            // check the SNP patterns for inconsistency
            checkBasePatternsForConsistency(potentialBasePatterns.keySet(), baseList, false);
        } while (baseListChanged);

        return baseList;
    }

    private static void updatePotentialBases(TreeSet basePattern, double proportion,
            TreeMap<Integer, HashMap> potentialBases, HashMap<Integer, SNP> baseList) {
        Iterator<SNP> SNP = basePattern.iterator();
        while (SNP.hasNext()) {
            SNP snp = SNP.next();

            if (baseList.containsKey(snp.getPosition())) { // position already assigned
                continue;
            }

            // get base proportions for this snp, if present
            HashMap<SNP, Double> baseProportions = potentialBases.get(snp.getPosition());
            if (baseProportions == null) { // this position has not been seen before
                // create a new map and add the proportion against this base
                baseProportions = new HashMap();
                baseProportions.put(snp, proportion);
                potentialBases.put(snp.getPosition(), baseProportions);

            } else {
                Double baseProportion = baseProportions.get(snp);
                // this base has not been seen at this position or better proportion found for this base
                if (baseProportion == null || proportion > baseProportion) {
                    baseProportions.put(snp, proportion);
                    potentialBases.put(snp.getPosition(), baseProportions);
                }
            }
        } // end while it.hasNext()

    }

    private static boolean assignBases(TreeMap<Integer, HashMap> potentialBases, LinkedHashMap<Integer, SNP> baseList) {

        boolean baseListChanged = false;
        Iterator<Map.Entry<Integer, HashMap>> itPotentialBases = potentialBases.entrySet().iterator();
        while (itPotentialBases.hasNext()) {
            // get the next entry
            Map.Entry<Integer, HashMap> entry = itPotentialBases.next();

            // get the position
            int position = entry.getKey();

            if (baseList.containsKey(position)) { // position already assinged
                continue;
            }

            // and the base proportions
            HashMap baseProportions = entry.getValue();
            Iterator<Map.Entry<SNP, Double>> itBaseProportions = baseProportions.entrySet().iterator();

            if (baseProportions.size() == 1) { // only one base found at this position
                // add it to the list
                baseList.put(position, itBaseProportions.next().getKey());
                // mark the base list as changed
                baseListChanged = true;
            } else {
                HashSet<SNP> selectedSNPs = new HashSet();
                double maximumProportion = 0.0;
                // assign the base with maximum  proprotion
                while (itBaseProportions.hasNext()) {
                    Map.Entry<SNP, Double> entryBP = itBaseProportions.next();

                    SNP snp = entryBP.getKey();
                    double proportion = entryBP.getValue();

                    if (proportion > maximumProportion) {
                        // save this snp as the selected snp;
                        selectedSNPs.clear();
                        selectedSNPs.add(snp);
                        // also update the maximum proportion
                        maximumProportion = proportion;
                    } else if (proportion == maximumProportion) {
                        selectedSNPs.add(snp);
                    }
                } // end while (itBaseProportions.hasNext())
                // add it to the list

                if (selectedSNPs.size() == 1) {
                    baseList.put(position, selectedSNPs.iterator().next());
                    // mark the base list as changed
                    baseListChanged = true;
                }
            }
        } // end while itPotentialBases.hasNext()

        return baseListChanged;
    }

    private static void checkBasePatternsForConsistency(Set<TreeSet> basePatterns, HashMap<Integer, SNP> baseList, boolean checkReverse) {

        if (baseList.isEmpty()) {
            return;
        }

        HashSet basePatternsToBeDeleted = new HashSet();
        // check all snps in this pattern
        Iterator<TreeSet> itSNPPattern = basePatterns.iterator();
        while (itSNPPattern.hasNext()) {
            // get the next SNP pattern
            TreeSet basePattern = itSNPPattern.next();

            // check if the base pattern is consistent with the already assigned bases
            boolean isConsistent = checkBasePattern(basePattern, baseList);
            if (checkReverse) {
                if (isConsistent) {
                    // Yes but we are checking for reverse. mark it for deletion
                    basePatternsToBeDeleted.add(basePattern);
                }
            } else {
                if (!isConsistent) {
                    // if not, mark it for deletion
                    basePatternsToBeDeleted.add(basePattern);
                }
            }
        } // end while it.hasNext()

        // remove the SNP patterns marked for deletion
        basePatterns.removeAll(basePatternsToBeDeleted);
    }

    private static boolean checkBasePattern(TreeSet basePattern, HashMap<Integer, SNP> baseList) {

        if (baseList.isEmpty()) {
            return true;
        }
        // check all snps in this pattern
        Iterator<SNP> itSNP = basePattern.iterator();
        while (itSNP.hasNext()) {
            // get the next SNP
            SNP theBase = itSNP.next();

            // get the SNP assigned at this position, if present
            SNP assignedBase = baseList.get(theBase.getPosition());
            if (assignedBase != null) { // only check assigned positions
                if (theBase.getBase() != assignedBase.getBase()) { // the base does not match the assigned base
                    return false;
                }
            }
        } // end while it.hasNext()

        return true;
    }

    private static HashSet<Integer> checkReferenceGenome(LinkedHashMap<Integer, SNP> geneHSPList,
            LinkedHashMap<Integer, SNP> geneReferenceBaseList, ArrayList<LinkedHashMap<Integer, SNP>> baseLists) {

        HashSet<Integer> referenceGenome = new HashSet<Integer>(DIPLOID_COUNT); // we use hashset as multiple genome might have the maximum identity with the reference

        // calculate identity of each genome with the reference
        double maxIdentity = 0.0;
        ArrayList<Double> identityList = new ArrayList<Double>(DIPLOID_COUNT);
        for (int d = 0; d < DIPLOID_COUNT; d++) {
            double identity = calculateReferenceIdentity(geneHSPList, geneReferenceBaseList, baseLists.get(d));
            identityList.set(d, identity);
            // also check for max identity
            maxIdentity = Math.max(identity, maxIdentity);
        }

        for (int d = 0; d < DIPLOID_COUNT; d++) {
            if (identityList.get(d) == maxIdentity) {
                referenceGenome.add(d);
            }
        }

        return referenceGenome;

    }

    private static double calculateReferenceIdentity(LinkedHashMap<Integer, SNP> geneHSPList,
            LinkedHashMap<Integer, SNP> geneReferenceBaseList, LinkedHashMap<Integer, SNP> baseList) {

        int identity = 0;
        // go through each HSP position for this gene
        Iterator<Integer> itHSPPosition = geneHSPList.keySet().iterator();
        while (itHSPPosition.hasNext()) {
            // get the next position
            int position = itHSPPosition.next();
            // get the base from the base list (if any) at the current position
            SNP theBase = baseList.get(position);
            if (theBase == null) { // snp was not assigned for this position, move to the next one
                continue;
            }
            SNP theReferenceBase = geneReferenceBaseList.get(position);

            // check if the base is same
            if (theReferenceBase.getBase() == theBase.getBase()) {
                identity++;
            }

        } // while (it.hasNext())

        return (double) identity / (double) geneHSPList.size();

    }

    private static void rectifyBaseListUsingReferenceGenome(LinkedHashMap<Integer, SNP> geneHSPList,
            LinkedHashMap<Integer, SNP> geneReferenceBaseList, LinkedHashMap<Integer, SNP> geneDiploidSNPList,
            LinkedHashMap<Integer, SNP> baseList) {

        // go through each HSP position for this gene
        Iterator<Integer> itHSPPosition = geneHSPList.keySet().iterator();
        while (itHSPPosition.hasNext()) {
            // get the next position
            int position = itHSPPosition.next();

            // get the base from the diploid base list at the current position
            SNP theDiploidBase = geneDiploidSNPList.get(position);
            if (theDiploidBase == null || !DNA.isNucleotide(theDiploidBase.getBase())) { // dont correct for positions where there is no coverage
                continue;
            }

            // get the base from the reference base list at the current position
            SNP theReferenceBase = geneReferenceBaseList.get(position);
            if (theReferenceBase.getBase() == DNA.BASE_N) { // ignore positions with N in the reference
                continue;
            }

            SNP theBase = baseList.get(position);
            if (theBase == null) {
                // snp was not assigned for this position or incorrectly assigned.
                // Replace it with reference base
                baseList.put(position, theReferenceBase);
            }
        } // while (it.hasNext())

    }

    private static void finaliseBaseList(LinkedHashMap<TreeSet, Double> basePatterns, HashMap<Integer, SNP> baseList, LinkedHashMap<TreeSet, Double> assignedSNPPatterns) {

        // assign genome SNPs
        HashMap newSNPList = assignBasesToGenome(basePatterns);

        // update the previous SNP list using the SNP list
        Iterator itSNPList = newSNPList.entrySet().iterator();
        while (itSNPList.hasNext()) {
            // get the next entry
            Map.Entry entry = (Map.Entry) itSNPList.next();
            // get the position
            int position = (Integer) entry.getKey();
            // and the SNP
            SNP newSNP = (SNP) entry.getValue();

            // get the SNP in the original SNP List
            SNP snp = (SNP) baseList.get(position);

            if (snp == null) { // snp not previously assigned
                baseList.put(position, newSNP);
                //AMT    
                //System.out.println("> Unassigned patterns used");
            } else if (!snp.equals(newSNP)) { // snp doesnt match the previously assigned SNP
                double maxIdentityAssigned = getMaximumIdentifyForBasePatterns(assignedSNPPatterns, snp);
                double maxIdentityNew = getMaximumIdentifyForBasePatterns(basePatterns, newSNP);
                if (maxIdentityNew > maxIdentityAssigned) {
                    baseList.put(position, newSNP);
                }

            }
        } // end while
    }

    private static double getMaximumIdentifyForBasePatterns(LinkedHashMap<TreeSet, Double> assignedPatterns, SNP theSNP) {

        double maximumIdentity = 0.0;
        // check all snps in this pattern
        Iterator<Map.Entry<TreeSet, Double>> itSNPPatterns = assignedPatterns.entrySet().iterator();
        while (itSNPPatterns.hasNext()) {
            Map.Entry<TreeSet, Double> entry = itSNPPatterns.next();
            // get the next SNP pattern
            TreeSet basePattern = entry.getKey();

            // check all snps in this pattern
            Iterator<SNP> itBasePattern = basePattern.iterator();
            while (itBasePattern.hasNext()) {
                // get the next SNP
                SNP currentSNP = itBasePattern.next();

                if (theSNP.equals(currentSNP)) {
                    maximumIdentity = Math.max(maximumIdentity, entry.getValue());
                    break;
                }
            } // end while itSNPPattern.hasNext()
        } // end while itSNPPatterns.hasNext()

        return maximumIdentity;

    }

    private static LinkedHashMap<Integer, Character> checkIfAllHSPBasesAreCovered(LinkedHashMap<Integer, SNP> geneHSPList, ArrayList<LinkedHashMap<Integer, SNP>> baseLists) {

        LinkedHashMap<Integer, Character> hspBasesCovered = new LinkedHashMap<Integer, Character>();

        // get all the bases assigned to sub-genomes at each position
        LinkedHashMap<Integer, TreeSet<Character>> consensusBaseList = new LinkedHashMap<Integer, TreeSet<Character>>();
        Iterator<LinkedHashMap<Integer, SNP>> itBaseLists = baseLists.iterator();
        while (itBaseLists.hasNext()) {

            LinkedHashMap<Integer, SNP> baseList = itBaseLists.next();
            Iterator<Integer> itPosition = baseList.keySet().iterator();
            while (itPosition.hasNext()) {
                int position = itPosition.next();

                TreeSet bases = consensusBaseList.get(position);
                if (bases == null) {
                    bases = new TreeSet();
                    consensusBaseList.put(position, bases);
                }
                bases.add(baseList.get(position).getBase());
            } // while (it.hasNext())
        }// while (itBaseLists.hasNext())

        // check if the consensus base from all assigned bases is same as HSP base
        Iterator<Integer> itPosition = geneHSPList.keySet().iterator();
        while (itPosition.hasNext()) {
            int position = itPosition.next();

            // get the HSP at this position
            SNP theHSP = (SNP) geneHSPList.get(position);

            // get the consensus base at this position
            TreeSet bases = consensusBaseList.get(position);
            if (bases == null) { //no base assigned at this position
                hspBasesCovered.put(position, 'N');
            } else {
                char consensusBase = (Character) DNA.getReverseNucleotideMap().get(bases.toString());
                if (theHSP.getBase() != consensusBase) {
                    hspBasesCovered.put(position, 'N');
                } else {
                    hspBasesCovered.put(position, 'Y');
                }
            }

        } // while (it.hasNext())
        return hspBasesCovered;
    }

    private static void writeBases(LinkedHashMap<Integer, SNP> baseList, LinkedHashMap<Integer, SNP> geneHSPList,
            LinkedHashMap<Integer, SNP> geneReferenceBaseList, LinkedHashMap<Integer, SNP> geneDiploidSNPList,
            String sequenceName, BufferedWriter out, LinkedHashMap<Integer, Character> hspBasesCovered) {

        Iterator<Integer> itPosition = geneHSPList.keySet().iterator();
        while (itPosition.hasNext()) {
            try {
                int position = itPosition.next();

                SNP theSNP = (SNP) baseList.get(position);
                if (theSNP == null) { // snp was not assigned for this genome, create a dummy snp
                    theSNP = new SNP(position, DNA.BASE_BLANK);
                }
                SNP theReferenceBase = geneReferenceBaseList.get(position);
                SNP theDiploidSNP = null;
                if (geneDiploidSNPList != null) {
                    theDiploidSNP = geneDiploidSNPList.get(position);

                    if (theDiploidSNP == null) {
                        theDiploidSNP = theReferenceBase;
                    } else if (theDiploidSNP.getBase() == DNA.BASE_BLANK) {
                        theDiploidSNP.setBase(theReferenceBase.getBase());
                    }
                } else {
                    theDiploidSNP = new SNP(position, DNA.BASE_AMBIGUOUS);
                }

                out.write(sequenceName + "\t" + position + "\t" + theReferenceBase.getBase() + "\t" + theDiploidSNP.getBase() + "\t" + theSNP.getBase() + "\t" + hspBasesCovered.get(position));
                out.newLine();

            } catch (IOException ex) {
                // do noting;
            }

        } // while (it.hasNext())
    }

    private static void printHelp() {

        System.out.println("HANDS v1.1");
        System.out.println("Usage: java -jar hands.jar <input parameters>");
        System.out.println("Input Parameters");
        System.out.println("\t-h or -help  :   Display this help");
        System.out.println("\t-s <str>     :   Polyploid SAM file");
        System.out.println("\t-g <str>     :   GFF file containing gene coordinates");
        System.out.println("\t-hsp <str>   :   Polyploid HSP file");
        System.out.println("\t-snp<n> <str>:   Diploid # n SNP file");
        System.out.println("\t-bd <str>    :   Polyploid Base Distribution File (Optional)");
        System.out.println("\t-bd<1> <str> :   Diploid # n Base Distribution File (Optional)");
        System.out.println("\t-out<n> <str>:   Sub-Genome # n Output File");
        System.out.println("\t-sp <double> :   SNP Pair Proportion Threshold (Default: 0.05)");
        System.out.println("\t-pm <double> :   SNP Pattern Matching Threshold (Default: 0.5)");
        System.out.println("\t-r <boolean> :   Rectify Assignment using Reference Genome (Default: FALSE)");
        System.out.println("\t-d <int>     :   Use Genome <int> as Distant Genome (Default: <null>)");
        System.out.println("Note: At most one Diploid SNP file can be missing. Use \"\" for the missing file.");
        System.out.println("      HANDS supports up to 10 genomes.");
    }
}
