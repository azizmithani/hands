/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package utility;

/**
 *
 * @author Aziz Mithani <aziz.mithani@lums.edu.pk>
 */
public class SAMAlignment {

    private String readName;
    private int flag;
    private String referenceSequenceName;
    private int alignmentStart;
    private int alignmentEnd;
    private int mappingQuality;
    private String cigar;
    private String mateReferenceSequenceName;
    private int mateStart;
    private int insertSize;
    private String readSequence;
    private String baseQualities;
    private String mismatchingPositions = "";
//    private String optionalFieldNM; // Edit distance
//    private String optionalFieldMD; // Mismatching positions/bases
//    private int optionalFieldAS; // Alignment score
//    private String optionalFieldBC; // Barcode sequence
//    private int optionalFieldX0; // Number of best hits
//    private int optionalFieldX1; // Number of suboptimal hits found by BWA
//    private int optionalFieldXN; // Number of ambiguous bases in the referenece
//    private int optionalFieldXM; // Number of mismatches in the alignment
//    private int optionalFieldXO; // Number of gap opens
//    private int optionalFieldXG; // Number of gap extentions
//    private String optionalFieldXT; // Type: Unique/Repeat/N/Mate-sw
//    private String optionalFieldXA; // Alternative hits; format: (chr,pos,CIGAR,NM;)*
//    private int optionalFieldXS; // Suboptimal alignment score
//    private String optionalFieldXF; // Support from forward/reverse alignment
//    private int optionalFieldXE; // Number of supporting seeds
    private SAMAlignment mate = null;

    
    /**
     * @return the readName
     */
    public String getReadName() {
        return readName;
    }

    /**
     * @param readName the readName to set
     */
    public void setReadName(String readName) {
        this.readName = readName;
    }

    /**
     * @return the flag
     */
    public int getFlag() {
        return flag;
    }

    /**
     * @param flag the flag to set
     */
    public void setFlag(int flag) {
        this.flag = flag;
    }


    /**
     * @return the alignmentStart
     */
    public int getAlignmentStart() {
        return alignmentStart;
    }

    /**
     * @param alignmentStart the alignmentStart to set
     */
    public void setAlignmentStart(int alignmentStart) {
        this.alignmentStart = alignmentStart;
    }

    /**
     * @return the alignmentEnd
     */
    public int getAlignmentEnd() {
        return alignmentEnd;
    }

    /**
     * @param alignmentEnd the alignmentEnd to set
     */
    public void setAlignmentEnd(int alignmentEnd) {
        this.alignmentEnd = alignmentEnd;
    }

    /**
     * @return the mappingQuality
     */
    public int getMappingQuality() {
        return mappingQuality;
    }

    /**
     * @param mappingQuality the mappingQuality to set
     */
    public void setMappingQuality(int mappingQuality) {
        this.mappingQuality = mappingQuality;
    }

    /**
     * @return the cigar
     */
    public String getCigar() {
        return cigar;
    }

    /**
     * @param cigar the cigar to set
     */
    public void setCigar(String cigar) {
        this.cigar = cigar;
    }

    /**
     * @return the mateStart
     */
    public int getMateStart() {
        return mateStart;
    }

    /**
     * @param mateStart the mateStart to set
     */
    public void setMateStart(int mateStart) {
        this.mateStart = mateStart;
    }

    /**
     * @return the insertSize
     */
    public int getInsertSize() {
        return insertSize;
    }

    /**
     * @param insertSize the insertSize to set
     */
    public void setInsertSize(int insertSize) {
        this.insertSize = insertSize;
    }

    /**
     * @return the readSequence
     */
    public String getReadSequence() {
        return readSequence;
    }

    /**
     * @param readSequence the readSequence to set
     */
    public void setReadSequence(String readSequence) {
        this.readSequence = readSequence;
    }

    /**
     * @return the baseQualities
     */
    public String getBaseQualities() {
        return baseQualities;
    }

    /**
     * @param baseQualities the baseQualities to set
     */
    public void setBaseQualities(String baseQualities) {
        this.baseQualities = baseQualities;
    }

    /**
     * @return the mate
     */
    public SAMAlignment getMate() {
        return mate;
    }

    /**
     * @param mate the mate to set
     */
    public void setMate(SAMAlignment mate) {
        this.mate = mate;
    }

    /**
     * @return the referenceSequenceName
     */
    public String getReferenceSequenceName() {
        return referenceSequenceName;
    }

    /**
     * @param referenceSequenceName the referenceSequenceName to set
     */
    public void setReferenceSequenceName(String referenceSequenceName) {
        this.referenceSequenceName = referenceSequenceName;
    }

    /**
     * @return the mateReferenceSequenceName
     */
    public String getMateReferenceSequenceName() {
        return mateReferenceSequenceName;
    }

    /**
     * @param mateReferenceSequenceName the mateReferenceSequenceName to set
     */
    public void setMateReferenceSequenceName(String mateReferenceSequenceName) {
        this.mateReferenceSequenceName = mateReferenceSequenceName;
    }

    /**
     * @return the mismatchingPositions
     */
    public String getMismatchingPositions() {
        return mismatchingPositions;
    }

    /**
     * @param mismatchingPositions the mismatchingPositions to set
     */
    public void setMismatchingPositions(String mismatchingPositions) {
        this.mismatchingPositions = mismatchingPositions;
    }

 
}
