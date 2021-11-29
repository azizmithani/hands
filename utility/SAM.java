/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package utility;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.HashMap;
import java.util.Iterator;
import java.util.LinkedHashSet;
import java.util.TreeSet;

/**
 *
 * @author Aziz Mithani <aziz.mithani@lums.edu.pk>
 */
public class SAM {

    // constants
    private final static int BASE_QUALITY_THRESHOLD = 20;
    public final static int FLAG_READ_PAIRED = 1;
    public final static int FLAG_MATE_UNMAPPED = 8;
    private static final char FIELDS_SEPERATOR = '\t';
    private static final char OPTIONAL_FIELD_SEPERATOR = ':';
    // variables
    String samFile;
    BufferedReader brSAM;
    private HashMap<String, Integer> sequenceLengths;

    public SAM(String samFile) {
        this.samFile = samFile;
        this.sequenceLengths = new HashMap();

        if (samFile.isEmpty()) {
            return;
        }

        try {
            this.brSAM = new BufferedReader(new FileReader(samFile));
        } catch (IOException ex) {
            // do nothing
        }

        readSequenceLengths();
    }

    private void readSequenceLengths() {

        try {
            while (peek(brSAM) == '@') {
                String line = brSAM.readLine();
                if (line.substring(1, 3).equals("SQ")) {
                    // get sequence lengths from the header
                    String[] entries = line.split("\t");

                    String sequenceName = entries[1].substring(3, entries[1].length());
                    int sequenceLength = Integer.parseInt(entries[2].substring(3, entries[2].length()));

                    sequenceLengths.put(sequenceName, sequenceLength);

                }

            } // end while
        } catch (IOException ex) {
            // do nothing
        }
    }

    private char peek(BufferedReader br) throws IOException {
        br.mark(1);
        char firstChar = (char) br.read();
        br.reset();

        return firstChar;
    }

    public SAMAlignment readNextAlignment() {

        String line;
        try {
            if ((line = brSAM.readLine()) != null) {

                SAMAlignment theSAMAlignment = new SAMAlignment();

                int pos = 0, end;
                // read name 
                end = line.indexOf(FIELDS_SEPERATOR, pos);
                theSAMAlignment.setReadName(line.substring(pos, end));
                //flag
                pos = end + 1;
                end = line.indexOf(FIELDS_SEPERATOR, pos);
                theSAMAlignment.setFlag(Integer.parseInt(line.substring(pos, end)));
                // sequence name 
                pos = end + 1;
                end = line.indexOf(FIELDS_SEPERATOR, pos);
                theSAMAlignment.setReferenceSequenceName(line.substring(pos, end));
                // read start
                pos = end + 1;
                end = line.indexOf(FIELDS_SEPERATOR, pos);
                theSAMAlignment.setAlignmentStart(Integer.parseInt(line.substring(pos, end)));
                // mapping quality
                pos = end + 1;
                end = line.indexOf(FIELDS_SEPERATOR, pos);
                theSAMAlignment.setMappingQuality(Integer.parseInt(line.substring(pos, end)));
                // cigar
                pos = end + 1;
                end = line.indexOf(FIELDS_SEPERATOR, pos);
                theSAMAlignment.setCigar(line.substring(pos, end));
                // mate sequence name 
                pos = end + 1;
                end = line.indexOf(FIELDS_SEPERATOR, pos);
                theSAMAlignment.setMateReferenceSequenceName(line.substring(pos, end));
                // mate start
                pos = end + 1;
                end = line.indexOf(FIELDS_SEPERATOR, pos);
                theSAMAlignment.setMateStart(Integer.parseInt(line.substring(pos, end)));
                // mate start
                pos = end + 1;
                end = line.indexOf(FIELDS_SEPERATOR, pos);
                theSAMAlignment.setInsertSize(Integer.parseInt(line.substring(pos, end)));
                // read sequence 
                pos = end + 1;
                end = line.indexOf(FIELDS_SEPERATOR, pos);
                theSAMAlignment.setReadSequence(line.substring(pos, end));
                // base qualities
                pos = end + 1;
                end = line.indexOf(FIELDS_SEPERATOR, pos);
                theSAMAlignment.setBaseQualities(line.substring(pos, end));
                // optional fields
                pos = end + 1;

                setOptionalFields(theSAMAlignment, line, pos);

                // calculate alignment end
                theSAMAlignment.setAlignmentEnd(getAlignmentEnd(theSAMAlignment.getAlignmentStart(), theSAMAlignment.getCigar()));

                return theSAMAlignment;
            } else {
                return null;
            }
        } catch (IOException ex) {
            // do nothing
            return null;
        }

    }

    /**
     * @return the sequenceLengths
     */
    public HashMap getSequenceLengths() {
        return sequenceLengths;
    }

    /**
     * @param sequenceLengths the sequenceLengths to set
     */
    public void setSequenceLengths(HashMap sequenceLengths) {
        this.sequenceLengths = sequenceLengths;
    }

    private void setOptionalFields(SAMAlignment theSAMAlignment, String alignment, int startPos) {

        int pos = startPos, end;
        do {
            // extract field name
            end = alignment.indexOf(OPTIONAL_FIELD_SEPERATOR, pos);
            String fieldName = alignment.substring(pos, end);
            // ignore field type
            pos = end + 1;
            end = alignment.indexOf(OPTIONAL_FIELD_SEPERATOR, pos);
            // extract field value
            pos = end + 1;
            end = alignment.indexOf(FIELDS_SEPERATOR, pos);
            String fieldValue = (end > 0 ? alignment.substring(pos, end) : alignment.substring(pos));
            // set pos value for next set of optional fields
            pos = end + 1;
            
            //TO-DO: Add other optional fields
            if (fieldName.equals("MD")) {
                theSAMAlignment.setMismatchingPositions(fieldValue);
            }
        } while (end > 0);


    }
    
    public static int getAlignmentEnd(int start, String cigar) {

        int sum = 0;
        String length = "";
        for (int i = 0; i < cigar.length(); i++) {
            if (Character.isDigit(cigar.charAt(i))) {
                length += cigar.charAt(i);
            } else if (cigar.charAt(i) == 'S') {
                // ignore
                length = "";
            } else if (cigar.charAt(i) == 'D') {
                sum += Integer.parseInt(length);
                length = "";
            } else if (cigar.charAt(i) == 'I') {
                length = "";
            } else if (cigar.charAt(i) == 'M' || cigar.charAt(i) == 'N') {
                sum += Integer.parseInt(length);
                length = "";
            }
        } // end for

        return start + sum - 1;
    }
    
    public static TreeSet<SNP> getMismatchingBases(LinkedHashSet<SNP> readRegionSNPList, int readStart, String MD, String readSeq, String baseQualities, String cigar) {

        // returns null if bad quality snps are present. also ignores SNPs found in deletion
        TreeSet<SNP> readSNPList = new TreeSet();
        Iterator<SNP> itReadRegionSNPList = readRegionSNPList.iterator();
        while (itReadRegionSNPList.hasNext()) {
            // get the next SNP
            SNP snp = itReadRegionSNPList.next();

            // Get the relative position of this SNP on the read
            int snpRelativePosition = getRelativePosition(snp.getPosition(), readStart, MD, cigar);

            // ignore this base if there's is a deletion in this read at SNP position (SNPRelativePosition == -1)
            if (snpRelativePosition == -1) {
                continue;
            }

            // return null this base if the sequencing quality is low
            int baseQuality = (int) (baseQualities.charAt(snpRelativePosition)) - 33; // base qualities are ascii - 33
            if (baseQuality < BASE_QUALITY_THRESHOLD) {
                return null;
            }

            // Get the base on the read
            char readBase = readSeq.charAt(snpRelativePosition);

            TreeSet bases = (TreeSet) DNA.getNucleotideMap().get(snp.getBase());
            // ignore this base if it's not part of the consensus
            if (bases == null || !bases.contains(readBase)) {
                continue;
            }

            // Add the SNP to this read's SNP list
            //readSNPList.add(Integer.toString(snpPosition) + ":" + readBase);
            readSNPList.add(new SNP(snp.getPosition(), readBase));

        }
        return readSNPList;
    }

    public static int getRelativePosition(int position, int readStart, String MD, String cigar) {

        // get the softclip present at the start of alignment, if any
        int softclip_length = 0;
        String softclip_length_str = "";
        for (int i = 0; i < cigar.length(); i++) {
            if (Character.isDigit(cigar.charAt(i))) {
                softclip_length_str += cigar.charAt(i);
            } else if (cigar.charAt(i) == 'S') {
                softclip_length += Integer.parseInt(softclip_length_str);
                break;
            } else if (cigar.charAt(i) == 'M' || cigar.charAt(i) == 'D'  || cigar.charAt(i) == 'I' || cigar.charAt(i) == 'N') {
                break;
            }
        } // end for

        // now parse the mismatch tag
        boolean inDeletion = false;
        String length = "";
        long sum = readStart;
        int relativePosition = position - readStart + softclip_length;
        int gap = 0;
        for (int i = 0; i < MD.length(); i++) {
            if (Character.isDigit(MD.charAt(i))) {
                if (inDeletion) {
                    inDeletion = false;
                }
                length += MD.charAt(i);
            } else if (MD.charAt(i) == 'A' || MD.charAt(i) == 'C' || MD.charAt(i) == 'G' || MD.charAt(i) == 'T') {
                if (!inDeletion) {
                    sum += Integer.parseInt(length);
                    if (sum >= position) {
                        break;
                    }
                } else {
                    if (sum == position) {
                        return -1;
                    }
                    gap++;
                }

                sum++; // increment the sum by one to include the current position
                length = "";
            } else if (MD.charAt(i) == '^') {
                sum += Integer.parseInt(length);
                if (sum > position) {
                    break;
                } else if (sum == position) {
                    return -1;
                }
                //RelativePosition += atoi(length.c_str());
                length = "";
                inDeletion = true;
            }

        } // end for
        return relativePosition - gap;
    }       
}
