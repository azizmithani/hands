/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package utility;

import java.util.HashMap;
import java.util.Iterator;
import java.util.Map;
import java.util.TreeSet;

/**
 *
 * @author Aziz Mithani <aziz.mithani@lums.edu.pk>
 */
public class DNA {

    public static final char BASE_BLANK = ' ';
    public static final char BASE_ZERO_COVERAGE = '0';
    public static final char BASE_LOW_COVERAGE = '<';
    public static final char BASE_AMBIGUOUS = '*';
    public static final char BASE_N = 'N';
    private static HashMap nucleotideMap;
    private static HashMap reverseNucleotideMap;

    private static void initialiseNucleotideMap() {

        TreeSet A = new TreeSet();
        A.add('A');
        TreeSet C = new TreeSet();
        C.add('C');
        TreeSet G = new TreeSet();
        G.add('G');
        TreeSet T = new TreeSet();
        T.add('T');
        TreeSet R = new TreeSet();
        R.add('A');
        R.add('G');
        TreeSet Y = new TreeSet();
        Y.add('C');
        Y.add('T');
        TreeSet M = new TreeSet();
        M.add('A');
        M.add('C');
        TreeSet K = new TreeSet();
        K.add('G');
        K.add('T');
        TreeSet S = new TreeSet();
        S.add('C');
        S.add('G');
        TreeSet W = new TreeSet();
        W.add('A');
        W.add('T');
        TreeSet H = new TreeSet();
        H.add('A');
        H.add('C');
        H.add('T');
        TreeSet B = new TreeSet();
        B.add('C');
        B.add('G');
        B.add('T');
        TreeSet D = new TreeSet();
        D.add('A');
        D.add('G');
        D.add('T');
        TreeSet V = new TreeSet();
        V.add('A');
        V.add('C');
        V.add('G');
        TreeSet N = new TreeSet();
        N.add('N');

        if (nucleotideMap == null) {
            nucleotideMap = new HashMap();
        }


        nucleotideMap.put('A', A);
        nucleotideMap.put('C', C);
        nucleotideMap.put('G', G);
        nucleotideMap.put('T', T);
        nucleotideMap.put('R', R);
        nucleotideMap.put('Y', Y);
        nucleotideMap.put('M', M);
        nucleotideMap.put('K', K);
        nucleotideMap.put('S', S);
        nucleotideMap.put('W', W);
        nucleotideMap.put('H', H);
        nucleotideMap.put('B', B);
        nucleotideMap.put('D', D);
        nucleotideMap.put('V', V);
        nucleotideMap.put('N', N);

    }

    private static void initialiseReverseNucleotideMap() {

        initialiseNucleotideMap();

        if (reverseNucleotideMap == null) {
            reverseNucleotideMap = new HashMap();
        }

        Iterator it = nucleotideMap.entrySet().iterator();
        while (it.hasNext()) {
            Map.Entry entry = (Map.Entry) it.next();
            // add the reverse pair
            reverseNucleotideMap.put(entry.getValue().toString(), entry.getKey());
        }

        // also add all four bases as 'N'
        TreeSet N = new TreeSet();
        N.add('A');
        N.add('C');
        N.add('G');
        N.add('T');
        reverseNucleotideMap.put(N.toString(), 'N');

    }

    public static HashMap getNucleotideMap() {

        if (nucleotideMap == null || nucleotideMap.isEmpty()) {
            initialiseNucleotideMap();
        }

        return nucleotideMap;
    }

    public static HashMap getReverseNucleotideMap() {

        if (reverseNucleotideMap == null || reverseNucleotideMap.isEmpty()) {
            initialiseReverseNucleotideMap();
        }

        return reverseNucleotideMap;
    }

    public static boolean isNucleotide(char base) {
        if (base == 'A' || base == 'T' || base == 'G' || base == 'C'
                || base == 'a' || base == 't' || base == 'g' || base == 'c') {
            return true;
        } else {
            return false;
        }
    }
    
    public static boolean isExtendedNucleotide(char base) {
    
        char upperBase = Character.toUpperCase(base);
        if (upperBase == 'R' || upperBase == 'Y' || upperBase == 'M' || upperBase == 'K' || upperBase == 'S' || upperBase == 'W' || upperBase == 'H' || upperBase == 'B' || upperBase == 'D' || upperBase == 'V' || upperBase == 'N') {
            return true;
        } else {
        return false;
    }
}
}
