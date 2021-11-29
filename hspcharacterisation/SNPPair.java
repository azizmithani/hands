/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package hspcharacterisation;

import utility.SNP;

/**
 *
 * @author Aziz Mithani <aziz.mithani@lums.edu.pk>
 */
public class SNPPair implements Comparable<Object> {

    private SNP snp1;
    private SNP snp2;
 
    @Override
    public int compareTo(Object object) {
        SNPPair snpPair = (SNPPair) object;

        if (this.snp1.compareTo(snpPair.snp1) < 0) {
            return -1;
        } else if (this.snp1.compareTo(snpPair.snp1) == 0) {
            if (this.snp2.compareTo(snpPair.snp2) < 0) {
                return -1;
            } else if (this.snp2.compareTo(snpPair.snp2) == 0) {
                return 0;
            } else {
                return 1;
            }
        } else {
            return 1;
        }
    }

    @Override
    public boolean equals(Object object) {
        if (!(object instanceof SNPPair)) {
            return false;
        }
        SNPPair snpPair = (SNPPair) object;
        return this.snp1.equals(snpPair.snp1) && this.snp2.equals(snpPair.snp2);
    }

    @Override
    public int hashCode() {
        int hash = 5;
        hash = 97 * hash + (this.snp1 != null ? this.snp1.hashCode() : 0);
        hash = 97 * hash + (this.snp2 != null ? this.snp2.hashCode() : 0);
        return hash;
    }

    public SNPPair(SNP snp1, SNP snp2) {//, int readCount) {
        this.snp1 = snp1;
        this.snp2 = snp2;
    //    this.readCount = readCount;
    }

    @Override
    public String toString() {
        return "(" + snp1 + "," + snp2 + ")";
    }

    public SNP getSNP1() {
        return snp1;
    }

    public SNP getSNP2() {
        return snp2;
    }


}