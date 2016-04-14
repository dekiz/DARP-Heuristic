/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package kapgelheuristic;

/**
 *
 * @author dileka
 */
public class Help {
    int[] bubble(float[] matrix) {
        int[] place = new int[matrix.length];
        for (int t = 0; t < place.length; t++) {
            place[t] = t;
        }
        for (int t1 = matrix.length - 1; t1 > 0; t1--) {
            for (int t2 = 0; t2 < t1; t2++) {
                if (matrix[t2] > matrix[t2 + 1]) {
                    float temp = matrix[t2];
                    int tplace = place[t2];
                    matrix[t2] = matrix[t2 + 1];
                    matrix[t2 + 1] = temp;
                    place[t2] = place[t2 + 1];
                    place[t2 + 1] = tplace;
                }
            }
        }
        return place;
    }//bubble ends
}
