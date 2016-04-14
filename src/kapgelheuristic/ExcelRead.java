/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package kapgelheuristic;

import java.io.File;
import jxl.Cell;
import jxl.Sheet;
import jxl.Workbook;
import org.apache.poi.ss.usermodel.DataFormatter;

/**
 *
 * @author dileka
 */
public class ExcelRead {
        public static String[] TimesRead(String excel) throws Exception {

//    File f=new File("C:\\Users\\IBM_ADMIN\\Documents\\NetBeansProjects\\KapgelProject\\KapgelProject\\simulated_data.xlsx");
        File f = new File(excel);
        Workbook wb = Workbook.getWorkbook(f);
        Sheet s = wb.getSheet(0);

        int col = s.getColumns();
        int row = s.getRows();

        double DestinationLat[] = new double[row];
        double DestinationLong[] = new double[row];
        double OriginLat[] = new double[row];
        double OriginLong[] = new double[row];
        String timestamp[] = new String[row];

        for (int j = 1; j < row; j++) {
            Cell A = s.getCell(2, j);
            String stringA = A.getContents();
//            System.out.println(stringA.substring(0, stringA.length()));
            DestinationLat[j] = Double.parseDouble(stringA);

            Cell B = s.getCell(3, j);
            String stringB = B.getContents();
            DestinationLong[j] = Double.parseDouble(stringB);
//            System.out.println(DestinationLong[j]);

            Cell C = s.getCell(4, j);
            String stringC = C.getContents();
//        System.out.println(stringA.substring(0, stringA.length()));
            OriginLat[j] = Double.parseDouble(stringC);

            Cell D = s.getCell(5, j);
            String stringD = D.getContents();
            OriginLong[j] = Double.parseDouble(stringD);

            // read timestamp
            Cell E = s.getCell(1, j);
            String stringE = E.getContents();
            timestamp[j] = stringE;
//            System.out.println(timestamp[j]);
        }

        String times[] = new String[row - 1];
        for (int i = 1; i < row; i++) {
            times[i - 1] = timestamp[i];
        }

        double coordinates[] = new double[4 * row - 4];
        for (int i = 1; i < row; i++) {
            coordinates[i - 1] = DestinationLat[i];
        }
        for (int i = 1; i < row; i++) {
            coordinates[i + row - 2] = DestinationLong[i];
        }
        for (int i = 1; i < row; i++) {
            coordinates[i + 2 * row - 3] = OriginLat[i];
        }

        for (int i = 1; i < row; i++) {
            coordinates[i + 3 * row - 4] = OriginLong[i];
        }

        for (int i = 0; i < coordinates.length; i++) {
//            System.out.println(coordinates[i]);
        }

        return times;

    }

     public static double[][] CoordinatesRead(String excel) throws Exception {
        DataFormatter formatter = new DataFormatter();
//    File f=new File("C:\\Users\\IBM_ADMIN\\Documents\\NetBeansProjects\\KapgelProject\\KapgelProject\\simulated_data.xlsx");

        File f = new File(excel);
        Workbook wb = Workbook.getWorkbook(f);
        Sheet s = wb.getSheet(0);

        int col = s.getColumns();
        int row = s.getRows();

        double DestinationLat[] = new double[row];
        double DestinationLong[] = new double[row];
        double OriginLat[] = new double[row];
        double OriginLong[] = new double[row];

        for (int j = 1; j < row; j++) {    
          
            Cell A = s.getCell(2,j);   
//            String stringA = A.getContents();
            String stringA = formatter.formatCellValue((org.apache.poi.ss.usermodel.Cell) A);
//            String tringA = A.getContents();
//            System.out.println(stringA.substring(0, stringA.length()));
//            DestinationLat[j] = Double.parseDouble(stringA);
            
            Cell B = s.getCell(3, j);
            String stringB = B.getContents();
            DestinationLong[j] = Double.parseDouble(stringB);
//            System.out.println(DestinationLong[j]);

            Cell C = s.getCell(4, j);
            String stringC = C.getContents();
//        System.out.println(stringC.substring(0, stringC.length()));
            OriginLat[j] = Double.parseDouble(stringC);

            Cell D = s.getCell(5, j);
            String stringD = D.getContents();
            OriginLong[j] = Double.parseDouble(stringD);

        }
        
        double pickupcoordinates[][] = {DestinationLat,DestinationLong};
        double [][] coordinates = new double[2*(pickupcoordinates[0].length-1)+2][2];
         coordinates[0][0] = 41.041999;
         coordinates[0][1] = 28.961532;
         coordinates[2*(pickupcoordinates[0].length-1)+1][0] = 41.041999;
         coordinates[2*(pickupcoordinates[0].length-1)+1][1] = 28.961532;
        for (int i = 1; i < pickupcoordinates[0].length; i++){
         for (int j=0; j<2; j++){
             coordinates[i][0] = DestinationLat[i];
             coordinates[i][1] = DestinationLong[i];
             coordinates[i+pickupcoordinates[0].length-1][0]=OriginLat[i];
             coordinates[i+pickupcoordinates[0].length-1][1]=OriginLong[i];

             
         }
     }
//        
        for (int i=0; i< (2*(pickupcoordinates[0].length-1)+2); i++){
            for (int j=0; j<2; j++){
                System.out.println(coordinates[i][j]);
            }
        }
        
        return coordinates;

    }
}
