/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package kapgelheuristic;

import java.util.*;
import java.io.*;

/**
 *
 * @author dileka
 */
public class KapgelHeuristic {

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {
        // TODO code application logic here
        long cputime = System.currentTimeMillis();
        long routetime = 0, rtime;
//Define parameters of the algorithm
        String file = "C:\\Users\\IBM_ADMIN\\Desktop\\Kapgel\\data\\d01"; //Name of data file
        final int M = 50; //Population size
        final int G = 15000; //Number of generations, iterations
        int Z = (int) M / 10; //Worst individuals in population
        if (Z == 0) {
            Z = 1; //At least 1
        }
        System.out.println("Z=" + Z + " M=" + M + " G=" + G + " " + file);
        final float pm = 1; //Mutation rate * 100
        final float big = 3.3E+38f; //big number
        float just = 60; //likely hood of choosing genes from p1
//float fjust=90; //max likely hood after iterations
//read data file + definations
        Request req = new Request(file);
        Car vehicle = new Car(file);
        int car = vehicle.getNOCar(); //number of cars available
        Route26042 ro = new Route26042(file);
        final int stops = req.getStops(); //number of stops
//I assume there are always at least 2 stops
//and there is an even number of them
        Chromo crowd = new Chromo(M, car, stops);
//make initial population
        crowd.createRandInitPop();
        byte[][][] population = crowd.getPopulation();
        byte[][][] population1 = new byte[M][car][stops / 2 + 1];
//Evaluate initial population
        byte[] clust = new byte[stops / 2 + 1]; //include all cust + depot
        float clustercost;
        float[] crowdcost = new float[M];
        float[] crowdcost1 = new float[M];
//cost for each individual in population
//System.out.println("clustercost ");
        for (int i = 0; i < M; i++) {
            clustercost = 0;
            for (int j = 0; j < car; j++) {
                clust = crowd.getCluster(i, j);
                rtime = System.currentTimeMillis();
                ro.getCustomers(clust, stops / 2);
                clustercost = clustercost + ro.getFinalcost();
                rtime = System.currentTimeMillis() - rtime;
                routetime = routetime + rtime;
            }//for ends
            crowdcost[i] = clustercost;
//System.out.print(crowdcost[i]+" ");
        }//for ends
//Genetic algorithm, G iterations
        byte[][] offspring = new byte[car][stops / 2 + 1];
        float offcost;
        float mutecost = big;
        for (int i = 0; i < G; i++) {
            if (i == 6000) {
                for (int tel1 = 0; tel1 < M; tel1++) {
                    crowdcost1[tel1] = crowdcost[tel1];
                    for (int s = 1; s < stops / 2 + 1; s++) {
                        int sum = 0, t = -1;
                        int[] ones = new int[car];
                        for (int c = 0; c < car; c++) {
                            population1[tel1][c][s] = population[tel1][c][s];
                            population1[tel1][c][0] = 1;
                            if (population[tel1][c][s] == 0) {
                                population[tel1][c][s] = 1;
                                t++;
                                ones[t] = c;
                            } else {
                                population[tel1][c][s] = 0;
                            }
                            sum = sum + population[tel1][c][s];
                        }//for
                        if (sum < 1) {
                            population[tel1][randInt(car)][s] = 1;
                        }
                        if (sum > 1) {
                            for (int p = 1; p < sum; p++) {
                                int rand = randInt(sum);
                                while (ones[rand] == -1) {
                                    rand = randInt(sum);
                                }
                                population[tel1][ones[rand]][s] = 0;
                                ones[rand] = -1;
                            }//for
                        }//if sum
                    }//for s
//Evaluate the new mutated population
                    clustercost = 0;
                    byte[] line = new byte[stops / 2 + 1];
                    for (int kk = 0; kk < car; kk++) {
                        for (int mm = 0; mm < stops / 2 + 1; mm++) {
                            line[mm] = population[tel1][kk][mm];
                        }
                        rtime = System.currentTimeMillis();
                        ro.getCustomers(line, stops / 2);
                        clustercost = clustercost + ro.getFinalcost();
                        rtime = System.currentTimeMillis() - rtime;
                        routetime = rtime + routetime;
                    }//for ends
                    crowdcost[tel1] = clustercost;
                }//for tel1
                printSolution(M, i, crowdcost);
            }//end (if i== 6000)
            if (i == 12000) {
                Help h = new Help();
                float[] cc1 = new float[crowdcost1.length];
                float[] cc2 = new float[crowdcost.length];
                for (int ix = 0; ix < crowdcost1.length; ix++) {
                    cc1[ix] = crowdcost1[ix];
                    cc2[ix] = crowdcost[ix];
                }
                int[] cost1 = h.bubble(cc1);
                int[] cost2 = h.bubble(cc2);
                for (int t2 = M / 2; t2 < cost2.length; t2++) {
                    int t1 = t2 - M / 2;
                    crowdcost[cost2[t2]] = crowdcost1[cost1[t1]];
                    population[cost2[t2]] = population1[cost1[t1]];
                }
            }//end if(i==12000)
//select individual stochastically from population
//as parent 1
            int rand;
            int parent1;
            rand = crowd.rand_int(100000000);
            parent1 = parenta(M, crowdcost, rand);
//select random individual from population as parent 2
//different than parent 1
            int parent2;
            do {
                parent2 = crowd.rand_int(M);
            } while (parent2 == parent1);
//if parent2 is better than parent1 then switch
            if (crowdcost[parent1] > crowdcost[parent2]) {
                int parentx;
                parentx = parent1;
                parent1 = parent2;
                parent2 = parentx;
            }
            byte[][] p1 = crowd.getIndividual(parent1);
            byte[][] p2 = crowd.getIndividual(parent2);
//offspring generated, crossover parent 1 and 2
            offspring
                    = crowd.cross_over2(p1, p2, just);
//evaluate offspring
            byte[] line = new byte[stops / 2 + 1];
            offcost = 0;
            for (int kk = 0; kk < car; kk++) {
                for (int mm = 0; mm < stops / 2 + 1; mm++) {
                    line[mm] = offspring[kk][mm];
                }
                rtime = System.currentTimeMillis();
                ro.getCustomers(line, stops / 2);
                offcost = offcost + ro.getFinalcost();
                rtime = System.currentTimeMillis() - rtime;
                routetime = rtime + routetime;
            }//for ends
//if offsprings cost same as parent 1 then mutate
            int key = 0;
            if (Math.abs(offcost - crowdcost[parent1]) < 1E-5f) {
                key = 80;
            }
//Mutation takes place with probability pm
            rand = crowd.rand_int(100);
            if (rand >= 0 && rand <= pm || key == 80) {
                key = 0;
                crowd.mutation(offspring);
//evaluate mutated offspring
                mutecost = 0;
                for (int kk = 0; kk < car; kk++) {
                    for (int mm = 0; mm < stops / 2 + 1; mm++) {
                        line[mm] = offspring[kk][mm];
                    }
                    rtime = System.currentTimeMillis();
                    ro.getCustomers(line, stops / 2);
                    mutecost = mutecost + ro.getFinalcost();
                    rtime = System.currentTimeMillis() - rtime;
                    routetime = rtime + routetime;
                }//for ends
            }//if ends
            if (mutecost < big) {
                offcost = mutecost;
            }
            mutecost = big;
//}//for j ends
//insert offspring for random individual among
//the worst in population
//Z is proportional to M
            float[] worst = new float[Z]; //cost of worst individuals
            int[] wno = new int[Z];
//indices of corresponding individuals
            for (int z = 0; z < worst.length; z++) {
                worst[z] = crowdcost[z];
                wno[z] = z;
            }//for ends
//small is the index of smallest element in worst
            int small = findSmall(worst, M);
            for (int t = Z; t < crowdcost.length; t++) {
                if (worst[small] < crowdcost[t]) {
                    worst[small] = crowdcost[t];
                    wno[small] = t;
                }
                small = findSmall(worst, M);
            }//for ends
//choose random individual among the worst with the higest cost
            int random = crowd.rand_int(wno.length);
            int b = wno[random];
//offspring replaces the chosen individual
//in the current population
            crowdcost[b] = offcost;
            for (int c = 0; c < car; c++) {
                for (int s = 0; s < stops / 2 + 1; s++) {
                    population[b][c][s] = offspring[c][s];
                }
            }
        }//for i ends
        cputime = System.currentTimeMillis() - cputime;
//print solution



        int small = findSmall(crowdcost, M);
        byte[] line = new byte[stops / 2 + 1];
        float[] ddist = new float[car];
        float[] samtalsridet = new float[car];
        float[] rodur = new float[car];
        float[] tw = new float[car];
        float[] waittime = new float[car];
        float[] bidtimifylki = new float[car];
        int customer = stops / 2;
        System.out.print(file + " " + customer);
        for (int kk = 0; kk < car; kk++) {
            for (int mm = 0; mm < stops / 2 + 1; mm++) {
                line[mm] = population[small][kk][mm];
            }
            ro.getCustomers(line, stops / 2);
            int sum = 0;
            int[] order = ro.getRoute();
            for (int k3 = 0; k3 < order.length; k3++) {
                sum = sum + order[k3];
                            for(int i=0; i<order.length; i++){
            System.out.println(order[i]);
            }
                            
            }
            int ld = 0;
            int[] load = ro.getLoad();
            for (int k6 = 0; k6 < load.length; k6++) {
                ld = ld + load[k6];
            }
            int leng = load.length;
            float prufa = leng - 2;
            System.out.print(" sum " + ld);
            System.out.print(" " + prufa);
            ddist[kk] = ro.getDist();
            float[] ride = ro.getRidetime();           
            for (int je = 0; je < ride.length; je++) {
                if (je % 2 == 1) {
                    samtalsridet[kk] = samtalsridet[kk] + ride[je];
                }
            }
            float[] time = ro.getTime();
            rodur[kk] = ro.getRouteDuration();
            float[] window = ro.getTimeWindowsViol();
            for (int k7 = 0; k7 < window.length; k7++) {
                tw[kk] = tw[kk] + window[k7];
            }
            float summa = 0f;
            float[] bidtimi = new float[load.length];
            float[] wait = ro.getWaitingTime();
            for (int k8 = 1; k8 < wait.length; k8++) {
                waittime[kk] = waittime[kk] + wait[k8];
                bidtimi[k8] = load[k8 - 1] * wait[k8];
                summa = summa + bidtimi[k8];
            }
            bidtimifylki[kk] = summa;
        }//for ends
        float samtals = 0;
        System.out.println("Distance for");
        for (int tel = 0; tel < ddist.length; tel++) {
System.out.print(" route "+tel+" is "+ddist[tel]);
            samtals = samtals + ddist[tel];
        }
                
        System.out.println();
        System.out.print(" dist " + samtals);
        System.out.println();
        samtals = 0;
        System.out.println("Ridetime for passangers in");
        for (int tel = 0; tel < samtalsridet.length; tel++) {
            System.out.print(" route " + tel + " is " + samtalsridet[tel]);
            samtals = samtals + samtalsridet[tel];
        }
        System.out.println();
        System.out.print(" ride " + samtals);
        System.out.println();
        samtals = 0;
        System.out.println("Route duration in");
        for (int tel = 0; tel < rodur.length; tel++) {
            System.out.print(" route " + tel + " is " + rodur[tel]);
            samtals = samtals + rodur[tel];
        }
        System.out.print(" route " + samtals);
        samtals = 0;
        System.out.println("Timewindows violation in");
        for (int tel = 0; tel < tw.length; tel++) {
            System.out.print(" route " + tel + " is " + tw[tel]);
            samtals = samtals + tw[tel];
        }
        System.out.println();
        System.out.print(" tw " + samtals);
        samtals = 0;
        System.out.println("Waiting times in");
        for (int tel = 0; tel < waittime.length; tel++) {
            System.out.print(" route " + tel + " is " + waittime[tel]);
            samtals = samtals + waittime[tel];
        }
        System.out.println();
        System.out.print(" wait " + samtals);
        System.out.println();
        samtals = 0;
        System.out.println("Waiting times in");
        for (int tel = 0; tel < bidtimifylki.length; tel++) {
            System.out.print(" route " + tel + " is " + waittime[tel]);
            samtals = samtals + bidtimifylki[tel];
        }
        System.out.println();
        System.out.print(" bid " + samtals);
        System.out.println();
        System.out.print(" M=" + M + " G=" + G + " pm=" + pm + " Z=" + Z + " " + file);
        System.out.print(" cost " + crowdcost[small]);
        System.out.print(" CPU " + cputime);
        System.out.println("min: " + cputime / 60000);
        System.out.print(" Routing running time " + routetime);
        System.out.println("Difference: " + (cputime - routetime) + " or ca "
                + ((cputime - routetime) * 100 / cputime) + "%");
        System.out.println("Number of customers is " + stops / 2);

        System.out.println();
        
    }

    //print solution
    //print solution
    public static void printSolution(int M, int i, float[] crowdcost) {
        int small = findSmall(crowdcost, M);
        System.out.print(i + " , " + crowdcost[small]);
        float[] cr = new float[crowdcost.length];
        for (int ij = 0; ij < cr.length; ij++) {
            cr[ij] = -crowdcost[ij];
        }
        small = findSmall(cr, M);
        System.out.println(" , " + crowdcost[small]);
        System.out.println();
        
    }

//finds smallest element in an array vec and returns
    public static int findSmall(int[] vec, int M) {
        int small = 2 * M, t;
        int sm = 2000000000, sm1;
        for (t = 0; t < vec.length; t++) {
            sm1 = sm;
            sm = Math.min(sm, vec[t]);
            if (sm != sm1) {
                small = t;
            }
        }//for ends
        return (small);
    }//findSmall ends

    public static int findSmall(float[] vec, int M) {
        int small = 2 * M, t;
        float sm = 3.0E+38f, sm1;
        for (t = 0; t < vec.length; t++) {
            sm1 = sm;
            sm = Math.min(sm, vec[t]);
            if (sm != sm1) {
                small = t;
            }
        }//for ends
        return (small);
    }//findSmall ends

    public static int parenta(int M, float[] cost, int rand) {
//select individual stochastically from population
//as parent 1
        float S = 0;
        int B = 100000000;
        for (int index = 0; index < M; index++) {
            S = S + 1 / cost[index];
        }
        int[] probability = new int[M];
        probability[0] = (int) (B / (cost[0] * S));
//probability*B of selecting ind 0 as parent 1
        int parent1 = M - 1; //default if no patent is chosen,
//because of numeration errors
        for (int index = 1; index < M; index++) {
            probability[index]
                    = probability[index - 1] + (int) (B / (cost[index] * S));
        }
        if (rand < probability[0]) {
            parent1 = 0;
        }
        for (int index = 1; index < M; index++) {
            if (rand >= probability[index - 1] && rand < probability[index]) {
                parent1 = index;
            }
        }//for ends
        return (parent1);
    } //parenta ends

    public static int parentb(int M, float[] cost, int rand) {
        int small = 0;
        float[] cro = new float[cost.length];
        for (int ij = 0; ij < cost.length; ij++) {
            cro[ij] = cost[ij];
        }
        for (int ij = 0; ij <= rand; ij++) {
            small = findSmall(cro, M);
            cro[small] = 3.3E36f;
        }
        return (small);
    }//parentb ends

    public static int randInt(int L) {
        return ((int) (Math.random() * L));
    }

    public static class Help {

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

    public static class Route26042 {

        private int[] v; //customers in cluster
        private int[] cs; //customers served
        private int[] cv; //customers in vehicle
        private int[] cn; //customers not jet served
        int MMM; //total number of customers
        String filename;
        int n; //total number of customers in all clusters
        int cun; //number of customers in cluster
        int speed = 1; //travelling time equal to E. dist
        float w1 = 1, w2; //(2)weight total route dur and time windows viol
        float w3, w4; //weight on ride time and route duration viol
        float w5 = 3, w6 = 1; //5(0)weight on ex ridetime and waiting time
        float w7 = 8; //weight on distance
        int[] ord; //array that holds the order of cust in route
        int[] load; //array with number of customers in vehicle after
//a node has been serviced
        Request r;
        Distance d;
        Car car;
        float[] ctime, wtime;

        public Route26042(String file) {
            filename = file;

        }

        public void giveOrder(int[] order, int nn, float begin) {
            ctime = new float[order.length];
            ctime[0] = begin;
            wtime = new float[order.length];
            load = new int[order.length];
            ord = new int[order.length];
            ord = order;
            cun = order.length / 2 - 1;
            n = nn;

        }
  
        public int[] getRoute() {
            return (ord);
            
        }

        public float getFinalcost() {
            if (ord.length > 2) {
                float finalcost = ordcost(ord);
                //System.out.println("finalcost " + finalcost);
                return (finalcost); 
                
            } else {
                return (0);
            }
        }

        public int[] getLoad() {
            return (load);
        }

        public float[] getTime() {
            if (ord.length > 2) {
                float cost = getFinalcost();
                System.out.println("Cost : "+cost);
            }
            return (ctime);
        }
//returns array with violated times

        public float[] getTimeWindowsViol() {
//float cost = getFinalcost();
            float[] out = new float[ord.length];
            for (int i = 0; i < ord.length; i++) {
                if (ctime[i] >= r.getLTimeWindow(ord[i]) + r.getServicetime(ord[i])
                        && ctime[i] <= r.getUTimeWindow(ord[i]) + r.getServicetime(ord[i])) {
                    out[i] = 0;
                } else if (ctime[i]
                        < r.getLTimeWindow(ord[i]) + r.getServicetime(ord[i])) {
                    out[i]
                            = ctime[i] - r.getLTimeWindow(ord[i]) - r.getServicetime(ord[i]);
                } else {
                    out[i]
                            = ctime[i] - r.getUTimeWindow(ord[i]) - r.getServicetime(ord[i]);
                }
                
            }//for ends
            return (out);
        }
//returns duration of a route

        public float getRouteDuration() {
            return (ctime[ctime.length - 1] - ctime[0]);
        }
//returns ride times of all customers in one route

        public float[] getRidetime() {
            float[] ride = new float[2 * cun];
            float cost = ordcost(ord);
            int tel = 0;
            for (int i1 = 1; i1 < ord.length - 1; i1++) {
                if (ord[i1] <= n) {
                    for (int i2 = 2; i2 < ord.length - 1; i2++) {
                        if (ord[i2] == ord[i1] + n) {
                            ride[tel] = ord[i1];
                            ride[tel + 1] = ctime[i2]
                                    - ctime[i1] - r.getServicetime(ord[i1]);
                            tel = tel + 2;
                        }//if
                    }//for
                }//if
            }//for
            return (ride);
        }//getRidetime
//returns waiting times of all nodes in one route

        public float[] getWaitingTime() {
            return wtime;
        }
//returns total distance between nodes in ord

        public float getDist() {
            float dist = 0;
            for (int i = 0; i < ord.length - 1; i++) {
                dist = dist + d.getDistance(ord[i], ord[i + 1]);
                }            
            return (dist);
            
        }
    
        
        
        public void getCustomers(byte[] clu, int MM) {
//clu inholds customers to be routed
            w2 = MM;
            w3 = MM;
            w4 = MM;
            cun = 0;
            for (int kk = 1; kk < clu.length; kk++) {
                cun = cun + clu[kk]; //number of customers in cluster
            }
            n = clu.length - 1; //total number of customers, depot not included
            r = new Request(filename);
            d = new Distance(filename);
            car = new Car(filename);
//if there are any customers in cluster then continue, otherwise
//return depot to depot route
            if (cun > 0) {
                v = new int[cun];
//v has the number of the customers in the cluster
                int count = 0;
                for (int i = 1; i < clu.length; i++) {
                    if (clu[i] == 1) {
                        v[count] = i;
                        count++;
                    }//if ends
                }//for ends
                route();
            }//if ends
            else {
                ord = new int[2];
                ord[0] = 0;
                ord[1] = 0;
                ctime = new float[2];
                ctime[0] = 0;
                ctime[1] = 0;
            }
        }//getCustomers ends
// -------------------- ROUTE ------------------------ //
//generate route for a cluster, algorithm taken from Baugh et al

        public void route() {
            ord = new int[2 * cun + 2]; //order of customers in solution
            load = new int[2 * cun + 2];
            ctime = new float[2 * cun + 2]; //time at which car has
//finished servicing the nodes,
            wtime = new float[2 * cun + 2];
            int nnode = -1, cnode = 0; //newnode and current node
            float ttime = 20000, mincost;
//order as in resulting route
            int sumcs = 0; //sum: cs - customers that have been serviced
            int sumcn = 0; //sum: cn - customers not in vehicle
            int sumv = 0; //sum: v - all customers numbers in cluster
            int[] firstcust = new int[v.length];
//first customer in route, cust no
            int[] firstno = new int[cun];
//number of customer in route
            cs = new int[cun];
            cv = new int[cun];
            cn = new int[cun];
            int[] N4 = new int[4]; //4 closest nodes to be considered as next node
//****************************************************
//find cust in v with earliest latest pick-up time
//as first customer in route
            int first = -5, no = -5;
            float mini = 3.3E38f;
            float[] ultw = new float[v.length];
            for (int re = 0; re < v.length; re++) {
                if (v[re] > n / 2) {
                    ultw[re] = r.getUTimeWindow(v[re]);
                } else {
                    ultw[re] = r.getUTimeWindow(v[re] + n)
                            - d.getDistance(v[re], v[re] + n) / speed
                            - r.getServicetime(v[re]);
                }
                if (mini > ultw[re]) {
                    mini = ultw[re];
                    no = re;
                    first = v[re];
                }
            }//for
//values set as -1 to indicate no customer for cs and cv, i.e.
//no customer is in vehicle or has been serviced in the begining,
//all customers number set into cn, customers not serviced
            for (int i = 0; i < v.length; i++) {
                cn[i] = v[i];
                cs[i] = -1;
                cv[i] = -1;
            }
            cv[no] = first; //first customer added to vehcile
            cn[no] = -1; //first customer deleted from cn
            ord[0] = 0; //start in depot
            ord[1] = first; //then to first customer
            load[0] = r.getLoadChange(0);
            load[1] = load[0] + r.getLoadChange(first);
//sums calulated, used as stopping criterias in while loop
            for (int i = 0; i < cs.length; i++) {
                sumcs = sumcs + cs[i];
                sumcn = sumcn + cn[i];
                sumv = sumv + v[i];
            }
//current node depot, next node first customer,
            cnode = 0;
            nnode = first;
            ttime = d.getDistance(cnode, nnode) / speed;
//travel time from depot to first customer
//find start time for depot, first cust has lowest upper pickup tw
            if (first > n / 2) {
                ctime[0] = r.getLTimeWindow(first) - ttime + r.getServicetime(0);
            } else {
                float L = r.getLTimeWindow(first + n)
                        - d.getDistance(first, first + n) / speed - r.getMaxRideTime();
                ctime[0] = L - ttime + r.getServicetime(0);
            }
            if (ctime[0] < 0) {
                ctime[0] = r.getServicetime(0);
            }
            ctime[1] = ctime[0] + ttime + r.getServicetime(nnode);
//time after first cust
//if ctime is lower than lower time window, then wait
            if (ctime[1] < r.getLTimeWindow(nnode) + r.getServicetime(nnode)) {
                ctime[1] = r.getLTimeWindow(nnode) + r.getServicetime(nnode);
                wtime[1] = ctime[1] - ctime[0] - ttime - r.getServicetime(nnode);
            }
//cnode set to first customer
            cnode = nnode;
            int kkk = 1; //counter in while loop, used in ord and ctime
//////////////// WHILE LOOP STARTS ///////////////////////
//for(int oj=0; oj<2*cun; oj++){//for perhaps better?
            while (sumcn > -cn.length || sumcs < sumv) {
                int cek = 0;
                for (int i = 0; i < N4.length; i++) {
                    N4[i] = -1;
                }
                for (int i = 0; i < N4.length; i++) {
//finds node closest to cnode but not included in N4
                    if (kkk < ctime.length) {
                        nnode = closest(cnode, ctime, N4, kkk);
                    } else {
                        System.out.println("L171 ctime problems");
                    }
//nnode=-2 when no more nodes are left
                    if (nnode == -2 && i == 0) {
                        System.out.println("KURT");
                    }
                    if (nnode == -2) {
                        break;
                    } else { //nnode set into N4
                        int count = 0;
                        for (int j = 0; j < N4.length; j++) {
                            if (N4[j] != -1) {
                                count++;
                            }
                        }
                        N4[count] = nnode;
                    }//else ends
                }//for ends
                mincost = 3.3E+38f;
//cost of nodes in N4 evaluated and cheapest chosen
                for (int i = 0; i < N4.length; i++) {
                    if (N4[i] > -1) {
                        if (mincost > cost(N4[i], cnode, ctime, kkk)) {
                            mincost = cost(N4[i], cnode, ctime, kkk);
                            nnode = N4[i];
                        }
                    }
                }
//Here something has gone wrong, algorithm should never have
//nnode= -2, it should stop before
                if (nnode == -2) {
                    System.out.println("L151 BREAK");
                    break;
                }//end if
                kkk++;
                visit(nnode);
                ttime = d.getDistance(cnode, nnode) / speed;
                if (kkk < ctime.length) {
                    ctime[kkk] = ctime[kkk - 1] + ttime + r.getServicetime(nnode);
                    if (ctime[kkk] < r.getLTimeWindow(nnode) + r.getServicetime(nnode)) {
                        ctime[kkk] = r.getLTimeWindow(nnode) + r.getServicetime(nnode);
                        wtime[kkk] = ctime[kkk] - r.getServicetime(nnode)
                                - ctime[kkk - 1] - ttime;
                    }
                }
                cnode = nnode;
                load[kkk] = load[kkk - 1] + r.getLoadChange(cnode);
//here the while loop has run more times than there are stops ??
                if (kkk > ord.length) {
                    System.out.println("NB!!!! EXTRA LOOP, cnode=" + cnode);
                    cek++;
                }
//if ok add cnode into route
                if (kkk < ord.length - 1 && cnode > 0) {
                    ord[kkk] = cnode;
                } else {
                    System.out.println("L289 Error in Route " + cnode);
                }
                if (cek > 0) {
                    System.out.println("Before: sumcs=" + sumcs + ", sumcn=" + sumcn);
                }
//recalculate stopping criterias
                sumcs = 0;
                sumcn = 0;
                for (int i = 0; i < cs.length; i++) {
                    sumcs = sumcs + cs[i];
                    sumcn = sumcn + cn[i];
                }
                if (cek > 0) {
                    System.out.println("After: sumcs=" + sumcs + ", sumcn=" + sumcn);
                    System.out.print(" sumv=" + sumv + " cn.length=" + cn.length);
                }
            }//while ends
        }//route ends
// ------------------- CLOSEST ----------------------- //
//Returns node closest to cnode

        public int closest(int cnode, float[] ctime, int[] N, int k) {
            float minldis = 3.3E+38f, newldis;
            int cstnode = -2, count;
//destination of c not in N
            for (int c = 0; c < cv.length; c++) {
                if (cv[c] != -1) {
                    count = 0;
                    for (int j = 0; j < N.length; j++) {
                        if (N[j] != cv[c] + n) {
                            count++;
                        }
                    }//for ends
                    if (count == N.length) {
                        newldis = separation(cnode, cv[c] + n, ctime, k);
                        if (minldis > newldis) {
                            minldis = newldis;
                            cstnode = cv[c] + n;
                        }//if ends
                    }//if ends
                }//if ends
            }//for ends
//now it is assumed that the customers travel alone,
//i.e. demand in each node is 1 or -1
            count = 0;
            for (int c = 0; c < cv.length; c++) {
                if (cv[c] != -1) {
                    count++;
                }
            }
            if (count == car.getCarCapacity()) {
                return (cstnode);
            }
            if (count > car.getCarCapacity()) {
                System.out.println("Capacity violated!!!!!!!");
                return (cstnode);
            }
//origin of c not in N
            for (int c = 0; c < cn.length; c++) {
                if (cn[c] != -1) {
                    count = 0;
                    for (int j = 0; j < N.length; j++) {
                        if (N[j] != cn[c]) {
                            count++;
                        }
                    }//for ends
                    if (count == N.length) {
                        newldis = separation(cnode, cn[c], ctime, k);
                        if (minldis > newldis) {
                            minldis = newldis;
                            cstnode = cn[c];
                        }//if ends
                    }//if ends
                }//if ends
            }//for ends
            return (cstnode);
        }//closest ends
// ---------------------- SEPARATION -------------------- //
//returns space-time separation between node cnode and nnode

        public float separation(int cnode, int nnode, float[] ctime, int k) {
            float ttime, timek1;
            float routedur = 0, twviol = 0, ridetimeviol = 0, exride = 0, wt = 0;
            ttime = d.getDistance(cnode, nnode) / speed;
            timek1 = ctime[k] + ttime; //arrival time at nnode
//change in route duration, ctime[k+1]-ctime[k]
            routedur = timek1 + r.getServicetime(nnode) - ctime[k];
//tw violations calculated
//for we want to get there soon (twviol>0)
            if (timek1 > r.getUTimeWindow(nnode)) {
                twviol = timek1 - r.getUTimeWindow(nnode);
            }
//tw viol - if early then increase cost - can wait(twviol<0)
            if (timek1 < r.getLTimeWindow(nnode)) {
                twviol = timek1 - r.getLTimeWindow(nnode);
            }
            if (timek1 < r.getLTimeWindow(nnode)) //if arrival is to early
            {
                timek1 = r.getLTimeWindow(nnode); //wait until ok
            }//customers ride time violations caluclated,
//if customer has been in the car for longer than max
//ride time sais then ride time viol > 0, service times of
//customer not included, i.e. nnode is a drop off location
// AND
//customers excess ride times calculated,
//ridetime - direct transport time
            if (nnode > n) {
                for (int i = 1; i < ord.length - 2; i++) {
                    if (ord[i] == nnode - n) {
                        if (timek1 - ctime[i] > r.getMaxRideTime()) {
                            ridetimeviol = timek1 - ctime[i];
                        }
                        exride
                                = timek1 - ctime[i] - d.getDistance(ord[i], nnode) / speed;
                    }//if
                }//for
            }//if
//waiting time calculated * persons in the vehicle
            int count = 0;
            for (int c = 0; c < cv.length; c++) {
                if (cv[c] != -1) {
                    count++;
                }
            }
            wt = count * (timek1 - (ctime[k] + ttime));
            return (w1 * routedur - w2 * twviol - w3 * ridetimeviol
                    - w5 * exride + w6 * wt + w7 * ttime * speed);
        }//separation ends
// ------------------------- COST ------------------------ //
//returns cost of visiting nnode

        public float cost(int nnode, int cnode, float[] ctime, int k) {
            int[] vn = new int[4]; //next four nodes considered
            int[] N = new int[4]; //empty array used in call to closest
            float totcost = 0;
            float ttime; //travel time
            int cc = 0;
            float ertime; //earliest arrival time to node
            for (int i = 0; i < 4; i++) {
                vn[i] = -1;
                N[i] = -1;
            }
            for (int i = 0; i < 4; i++) {
                totcost = totcost + movecost(cnode, nnode, ctime, k);
                if (nnode > 0) {
                    visit(nnode);
                }
                vn[i] = nnode;
                ttime = d.getDistance(cnode, nnode) / speed;
                ertime = ctime[k] + ttime; //arrival time at nnode
                if (ertime < r.getLTimeWindow(nnode)) {
                    ertime = r.getLTimeWindow(nnode);
                }
                cnode = nnode;
                k++;
                ctime[k] = ertime + r.getServicetime(cnode);
//service at nnode finished
                nnode = closest(cnode, ctime, N, k);
                if (nnode == -2) {
                    cc++;
                    if (cc == 1) {
                        nnode = 0;
                        totcost = totcost + movecost(cnode, nnode, ctime, k);
                    } else {
                        break;
                    }
                }//end if
            }//end for
            for (int i = vn.length - 1; i >= 0; i--) {
                if (vn[i] > 0) {
                    unvisit(vn[i]);
                }
            }//for ends
            return (totcost);
        }//cost ends
// ----------------------- MOVECOST ------------------- //
//returns cost of move from current node to next node

        public float movecost(int cnode, int nnode, float[] ctime, int k) {
            float ttime, timek1;
            float routedur = 0, twviol = 0, ridetimeviol = 0, exride = 0, wt = 0;
            ttime = d.getDistance(cnode, nnode) / speed;
            timek1 = ctime[k] + ttime; //arrival time at nnode
//change in route duration, ctime[k+1]-ctime[k]
            routedur = timek1 + wtime[k + 1] + r.getServicetime(nnode) - ctime[k];
//tw violations calculated,
            if (timek1 > r.getUTimeWindow(nnode)) {
                twviol = timek1 - r.getUTimeWindow(nnode);
            }
            if (ctime[k] + ttime < r.getLTimeWindow(nnode)) {
                twviol = r.getLTimeWindow(nnode) - (ctime[k] + ttime);
            }
            if (timek1 < r.getLTimeWindow(nnode)) //if arrival is to early
            {
                timek1 = r.getLTimeWindow(nnode); //wait until ok
            }//customers ride time violations caluclated,
//if customer has been in the car for longer than max
//ride time sais then ride time viol > 0, service times of
//customer not included, i.e. nnode is a drop off location
// AND
//customers excess ride times calculated,
//ridetime - direct transport time
            if (nnode > n) {
                for (int i = 1; i < ord.length - 2; i++) {
                    if (ord[i] == nnode - n) {
                        if (timek1 - ctime[i] > r.getMaxRideTime()) {
                            ridetimeviol = timek1 - ctime[i] - r.getMaxRideTime();
                        }
                        exride
                                = timek1 - ctime[i] - d.getDistance(ord[i], nnode) / speed;
                    }//if
                }//for
            }//if
//waiting time calculated * persons in the vehicle
            int count = 0;
            for (int c = 0; c < cv.length; c++) {
                if (cv[c] != -1) {
                    count++;
                }
            }
            wt = count * (timek1 - (ctime[k] + ttime));
            return (w1 * routedur + w2 * twviol + w3 * ridetimeviol
                    + w5 * exride + w6 * wt + w7 * ttime * speed);
        }//movecost ends
// ----------------------- ORDCOST ------------------- //
//Returns cost of route

        public float ordcost(int[] order) {
            float ttime, routeviol = 0, serv, routedur = 0, twviol = 0;
            float rideviol = 0, xride = 0, wt = 0;
            int cnode, nnode;
            for (int i = 0; i < cv.length; i++) {
                cv[i] = -1;
                cs[i] = -1;
                cn[i] = v[i];
            }
            load[0] = 0;
            cnode = order[0];
            ctime[0] = retime();
            for (int j = 0; j < wtime.length; j++) {
                wtime[j] = 0;
            }
            
//update ctime, wtime, load for the route found
            for (int i = 1; i < order.length; i++) {
                nnode = order[i];
                load[i] = load[i - 1] + r.getLoadChange(order[i]);
                ttime = d.getDistance(cnode, nnode) / speed;
                ctime[i] = ctime[i - 1] + ttime + r.getServicetime(nnode) + wtime[i];
                if (ctime[i] < r.getLTimeWindow(nnode) + r.getServicetime(nnode)) {
                    ctime[i] = r.getLTimeWindow(nnode) + r.getServicetime(nnode);
                    wtime[i] = ctime[i] - (ctime[i - 1] + ttime + r.getServicetime(nnode));
                }//if ends
                if (i > 1 && wtime[i] > 0 && load[i - 1] > load[i - 2]) {
                    wait(i);
                }
//if(diff>-1) totcost = diff;
                cnode = nnode;
                if (cnode > 0) {
                    visit(cnode);
                }
            }//for ends
            for (int i = 1; i < order.length; i++) {
                ttime = d.getDistance(order[i - 1], order[i]) / speed;
                serv = r.getServicetime(order[i]);
//tw violations calculated,
                if (ctime[i] - serv > r.getUTimeWindow(order[i])) {
                    twviol = twviol + (ctime[i] - serv) - r.getUTimeWindow(order[i]);
                }
                if (ctime[i - 1] + ttime < r.getLTimeWindow(order[i])) {
                    twviol = twviol + r.getLTimeWindow(order[i]) - (ctime[i - 1] + ttime);
                }
//customers ride time violations calculated,
//if customer has been in the car for longer than max
//ride time sais then ride time viol > 0, service times of
//customer not included, i.e. nnode is a drop off location
// AND
//customers excess ride times calculated,
//ridetime - direct transport time
                if (order[i] > n) {
                    for (int j = 1; j < order.length - 2; j++) {
                        if (order[j] == order[i] - n) {
                            if (ctime[i] - serv - ctime[j] > r.getMaxRideTime()) {
                                rideviol = rideviol + ctime[i] - serv
                                        - ctime[j] - r.getMaxRideTime();
                            }
                            xride = xride + ctime[i] - serv - ctime[j]
                                    - d.getDistance(order[i], order[j]) / speed;
                        }//if
                    }//for
                }//if
//waiting time calculated * persons in the vehicle
                wt = wt + load[i - 1] * (ctime[i] - serv - (ctime[i - 1] + ttime));
            }
            if (wt < 1E-5) {
                wt = 0;
            }
//total route duration
            routedur = ctime[ctime.length - 1] - ctime[0];
//route duration violations calculated,
//if current time is higer than
//max allowable route duration then route viol becomes > 0
            if (ctime[ctime.length - 1] - ctime[0] > car.getRouteDuration()) {
                routeviol
                        = ctime[ctime.length - 1] - ctime[0] - car.getRouteDuration();
            }
            return (w1 * routedur + w2 * twviol + w3 * rideviol
                    + w4 * routeviol + w5 * xride + w6 * wt + w7 * getDist());
        }//ordcost ends
        
// ------------------------- RETIME -------------------------//
//Recalculate starting time so that time windows are not violated

        public float retime() {
            float ttime;
            if (ctime[ctime.length - 1] > r.getUTimeWindow(ord[ord.length - 1])) {
                ctime[ctime.length - 1] = r.getUTimeWindow(ord[ord.length - 1]);
            }
            for (int i = ctime.length - 1; i > 0; i--) {
                ttime = d.getDistance(ord[i], ord[i - 1]) / speed;
                ctime[i - 1] = ctime[i] - ttime - r.getServicetime(ord[i]);
                if (ctime[i - 1]
                        < r.getLTimeWindow(ord[i - 1]) + r.getServicetime(ord[i - 1])) {
                    ctime[i - 1]
                            = r.getLTimeWindow(ord[i - 1]) + r.getServicetime(ord[i - 1]);
                }
                if (ctime[i - 1]
                        > r.getUTimeWindow(ord[i - 1]) + r.getServicetime(ord[i - 1])) {
                    ctime[i - 1]
                            = r.getUTimeWindow(ord[i - 1]) + r.getServicetime(ord[i - 1]);
                }
            }//for ends
            return (ctime[0]);
        }
        

// -------------------------- WAIT -----------------------//
//moves waiting times to nodes where there are fewer customers
//in the vehicle

        public void wait(int no) {
            float start = wtime[no];
            float diff = 0;
            for (int i = no; i > 1; i--) {
                diff = r.getUTimeWindow(ord[i - 1])
                        - ctime[i - 1] + r.getServicetime(ord[i - 1]);
//arrival time within tw
                if (load[i - 2] > load[i - 1] || wtime[i] < 0.01 || diff <= 0) {
                    break;
                } else {
                    if (diff > 0 && diff <= wtime[i]) {
                        wtime[i - 1] = wtime[i - 1] + diff;
                        ctime[i - 1] = ctime[i - 1] + diff;
                        wtime[i] = wtime[i] - diff;
                        diff = 0;
                    }
                    if (diff > 0 && diff > wtime[i]) {
                        wtime[i - 1] = wtime[i - 1] + wtime[i];
                        ctime[i - 1] = ctime[i - 1] + wtime[i];
                        wtime[i] = 0;
                    }
                }//else ends
            }//for ends
        }//wait ends
// -------------------------- VISIT -------------------- //
//marks gnode as visited by updating global data

        public void visit(int gnode) {
            for (int c = 0; c < v.length; c++) {
                if (gnode == v[c]) {
                    cn[c] = -1;
                    cv[c] = v[c];
                    break;
                }
                if (gnode == v[c] + n) {
                    cs[c] = v[c];
                    cv[c] = -1;
                    break;
                }
            }//for ends
        }//visit ends
// -------------------------- UNVISIT -------------------- //
//marks gnode as not visited by updating global data

        public void unvisit(int gnode) {
            for (int c = v.length - 1; c >= 0; c--) {
                if (gnode == v[c]) {
                    cn[c] = v[c];
                    cv[c] = -1;
                    break;
                }
                if (gnode == v[c] + n) {
                    cs[c] = -1;
                    cv[c] = v[c];
                    break;
                }
            }//for ends
        }//uvisit ends
    }//class ends

    public static class Request {

        int[][] req;
        float[][] coo;
        Data s;
        int st, max;

        public Request() {
        }
//Constructer that initializes matrices and reads in values

        public Request(String filename) {
            final int sizereq = 4, sizecoo = 2;
            Data s = new Data(filename);
            st = s.getStops();
            req = new int[st][sizereq];
            coo = new float[st][sizecoo];
            coo = s.getCoo();
            req = s.getReq();
            max = s.getMaxRideTime();
        }//Request constructer ends

        public int getStops() {
            return (st);
        }

        public int getMaxRideTime() {
            return (max);
        }

        public float[][] getCooMatrix() {
            return (coo);
        }//getCooMatrix ends

        public float getxCoordinates(int custnum) {
            float x = coo[custnum][0];
            return (x);
        }//getxCoordinates ends

        public float getyCoordinates(int custnum) {
            float y = coo[custnum][1];
            return (y);
        }//getyCoordinates ends

        public float[] getxy(int custnum) {
            float[] xy = new float[2];
            xy[0] = getxCoordinates(custnum);
            xy[1] = getyCoordinates(custnum);
            return (xy);
        }

        public int getServicetime(int custnum) {
            int s = req[custnum][0];
            return (s);
        }//getServicetime ends

        public int getLoadChange(int custnum) {
            int c = req[custnum][1];
            return (c);
        }//getLoadChange ends

        public int getLTimeWindow(int custnum) {
            int tw = req[custnum][2];
            return (tw);
        }//getlTimeWindow ends

        public int getUTimeWindow(int custnum) {
            int tu = req[custnum][3];
            return (tu);
        }//getuTimeWindow ends
    }//Request class ends

    public static class Data {

        private LineNumberReader in;
        int[] prob = new int[5];
        float[][] coo;
        int[][] req;
        int count;
        String filename;
//method that reads number of depots, number of stops, max duration
//time of a tour, max capacity of cars, max riding time of customers

        public Data(String file) {
            filename = file;
            readProblem();
            read();
        }

        public void readProblem() {
            try {
                in = new LineNumberReader(new FileReader(filename));
                StringTokenizer dimen;
//reads the information into a vector prob
                for (int i = 1; i < 2; i++) {
                    String dimension = in.readLine();
                    dimen = new StringTokenizer(dimension);
                    while (dimen.hasMoreTokens()) {
                        prob[0] = Integer.parseInt(dimen.nextToken());
                        prob[1] = Integer.parseInt(dimen.nextToken());
                        prob[2] = Integer.parseInt(dimen.nextToken());
                        prob[3] = Integer.parseInt(dimen.nextToken());
                        prob[4] = Integer.parseInt(dimen.nextToken());
                    }//while ends
                }//for ends
            }//try ends
            catch (EOFException eof) {
                closeFile();
            } catch (IOException e) {
                System.out.println("1 The file " + filename + " could not be opened "
                        + e.toString());
                System.exit(1);
            }//CATCH ENDS
        }//readProblem ends
//method that returnes number of vechicles available

        int getNOCar() {
            return (prob[0]);
        }
//method that returnes number of stops

        int getStops() {
            return (prob[1]);
        }
//method that returns allowable route duration

        int getRouteDuration() {
            return (prob[2]);
        }
//method that returns capacity of cars

        int getCapacity() {
            return (prob[3]);
        }
//method that returns maximum riding time for customers

        int getMaxRideTime() {
            return (prob[4]);
        }
//method that reads cooridnates of customers into a matrix coo and
//service time, load change, lower time window and upper time window
//into matrix req

        public void read() {
            try {
                in = new LineNumberReader(new FileReader(filename));
                StringTokenizer tokens;
//Reads cooridinates of the customers into a matrix coo
//and rest into a matrix req
                int dim = getStops() + 1;
                int car = getNOCar();
                req = new int[dim][4];
                coo = new float[dim][2];
                int d = 0;
                if (dim - 1 > 0 && dim - 1 < 10) {
                    d = 0;
                }
                if (dim - 1 > 9 && dim - 1 < 100) {
                    d = 1;
                }
                if (dim - 1 > 99 && dim - 1 < 1000) {
                    d = 2;
                }
                if (car > 9 && car < 100) {
                    d = d + 1;
                }
                in.skip(12 + d);
                for (int i = 0; i < dim + 1; i++) {
                    String tokenstring = in.readLine();
                    int index;
                    tokens = new StringTokenizer(tokenstring);
                    while (tokens.hasMoreTokens()) {
                        index = Integer.parseInt(tokens.nextToken());
                        coo[index][0] = Float.parseFloat(tokens.nextToken());
                        coo[index][1]
                                = Float.parseFloat(tokens.nextToken());
                        req[index][0]
                                = Integer.parseInt(tokens.nextToken());
                        req[index][1]
                                = Integer.parseInt(tokens.nextToken());
                        req[index][2]
                                = Integer.parseInt(tokens.nextToken());
                        req[index][3]
                                = Integer.parseInt(tokens.nextToken());
                    }//WHILE ENDS
                
                
                } //FOR ENDS
            } //TRY ENDS
            catch (EOFException eof) {
                closeFile();
            } catch (IOException e) {
                System.out.println("2 The file " + filename + " could not be opened "
                        + e.toString());
            }//CATCH ENDS
        }//readRequest ends
    
    
        
        private void closeFile() {
            try {
                in.close();
                System.exit(0);
            } catch (IOException e) {
                System.err.println("Error closing file" + e.toString());
                System.exit(1);
            }
        }//closeFile ends

        public float[][] getCoo() {
            return (coo);
                    }

        public int[][] getReq() {
            return (req);
        }
    } //CLASS Data ENDS

    public static class Distance {

        protected float[][] dist;
        int sto;
        float[][] dcoo;

        public Distance(String filename) {
            Request r = new Request(filename);
            sto = r.getStops();
            dist = new float[sto + 1][sto + 1];
            dcoo = r.getCooMatrix();
            calculateDistance();
        }//Distance Constructur ends
        
        public static final double R = 6372.8;
        
        public void calculateDistance() {
            
            for (int s = 0; s < sto + 1; s++) {
                dist[s][s] = 0; 
            }
            
            for (int s1 = 0; s1 < sto + 1; s1++) {
                for (int s2 = 0; s2 < s1; s2++) {
                    double dLat = Math.toRadians(dcoo[s2][1] - dcoo[s1][1]);
                    double dLon = Math.toRadians(dcoo[s2][0] - dcoo[s1][0]);
                    double lat1=Math.toRadians(dcoo[s1][1]);
                    double lat2=Math.toRadians(dcoo[s1][1]);
                    
                    double a = Math.pow(Math.sin(dLat / 2),2) + Math.pow(Math.sin(dLon / 2),2) * Math.cos(lat1) * Math.cos(lat2);
                    double c = R * 2 * Math.asin(Math.sqrt(a));
                           dist[s1][s2] = (float) c;              
//                    dist[s1][s2]
//                            = (float) Math.sqrt((dcoo[s1][0] - dcoo[s2][0])
//                                    * (dcoo[s1][0] - dcoo[s2][0])
//                                    + (dcoo[s1][1] - dcoo[s2][1])
//                                    * (dcoo[s1][1] - dcoo[s2][1]));
                   
                }//for ends
                
            }//for ends
            for (int s1 = 0; s1 < sto + 1; s1++) {
                for (int s2 = s1 + 1; s2 < sto + 1; s2++) {
                    dist[s1][s2] = dist[s2][s1];
                }//for ends
            }//for ends
        }//calculateDistance ends
//returnes distance between stop 1 and 2

        public float getDistance(int cust1, int cust2) {          
            return (dist[cust1][cust2]);
        }//getDistance ends
    }//Distance class ends

    public static class Car {

        int capacity, carno;
        String filename;
        Data d;

        public Car(String file) {
            filename = file;
            d = new Data(filename);
        }

        public int getCarCapacity() {
            return d.getCapacity();
        }

        public int getNOCar() {
            return d.getNOCar();
        }

        public int getRouteDuration() {
            return d.getRouteDuration();
        }
    }//Car class ends

    public static class Chromo extends Request {

        byte[][][] ind;
        int[][] sum;
        int stops, car, pop;

        public Chromo() {
        }

        public Chromo(int p, int c, int st) {
            pop = p;
            car = c;
            stops = st;
        }

        public void createRandInitPop() {
            int rand;
            ind = new byte[pop][car][stops / 2 + 1];
            for (int i = 0; i < pop; i++) {
                for (int j = 0; j < car; j++) {
                    ind[i][j][0] = 1;
                }
                for (int k = 1; k < stops / 2 + 1; k++) {
                    rand = rand_int(car);
                    ind[i][rand][k] = 1;
                }//for ends
            }//for ends
            checkPopulation();
//printMatrix(ind[0]);
//printMatrix(ind[1]);
        }//createRandIndPop ends

        public byte[][][] getPopulation() {
            return ind;
        }

        public byte[][] getIndividual(int a) {
            byte[][] parent = new byte[car][stops / 2 + 1];
            for (int i = 0; i < car; i++) {
                for (int j = 0; j < stops / 2 + 1; j++) {
                    parent[i][j] = ind[a][i][j];
                }
            }
            return (parent);
        }

        public byte[] getCluster(int a, int b) {
            byte[] clust = new byte[stops / 2 + 1];
            for (int i = 0; i < stops / 2 + 1; i++) {
                clust[i] = ind[a][b][i];
            }
            return (clust);
        }

        public void checkPopulation() {
            int rand;
            for (int i = 0; i < pop; i++) {
                int sum = 0;
                for (int j = 0; j < car; j++) {
                    sum = sum + ind[i][j][0];
                }
                if (sum != car) {
                    System.out.println("Error for depot " + sum + " car " + car);
                }
                for (int k = 1; k < stops / 2 + 1; k++) {
                    int su = 0;
                    for (int j = 0; j < car; j++) {
                        su = su + ind[i][j][k];
                    }
                    switch (su) {
                        case 0:
                            rand = rand_int(car);
                            ind[i][rand][k] = 1;
                            break;
                        case 1:
                            break;
                        default:
                            System.out.println("Error individiual " + i + " car " + k);
                    }//switch ends
                }//for ends
            }//for ends
        }//checkPopulation ends

        public void custInCar() {
            sum = new int[pop][car];
            for (int i = 0; i < pop; i++) {
                for (int j = 0; j < car; j++) {
                    for (int k = 0; k < stops / 2; k++) {
                        sum[i][j] = sum[i][j] + ind[i][j][k];
                    }//for ends
                }//for ends
            }//for ends
        }//CustInCar ends

        public int getCustInCar(int ind, int carno) {
            custInCar();
            return (sum[ind][carno]);
        }
//cross over in which there is selected randomly rows from
//both parents, random template created and cross over done
//accordingly. Offspring is then same as first parent but
//row chosen replaced by cross over line. Then the offspring
//is checked to see if all customers are included once, if
//not the error is corrected, tried to keep the cross over
//line unchanged, to avoid duplicates: offspring = parent 1

        public byte[][] cross_over1(byte[][] inda, byte[][] indb) {
            int row = inda.length; //# rows
            int col = inda[0].length; //# columns in row 0
            byte[][] offspring = new byte[row][col];
            for (int i = 0; i < row; i++) {
                for (int j = 0; j < col; j++) {
                    offspring[i][j] = inda[i][j];
                }
            }
            int randa = rand_int(row);
            int randb = rand_int(row);
            int[] template = new int[col];
            for (int i = 0; i < col; i++) {
                template[i] = rand_int(2);
            }
            for (int i = 0; i < col; i++) {
                if (template[i] == 1) {
                    offspring[randa][i] = indb[randb][i];
                }
            }//for ends
            offspring = correct_Matrix(offspring, randa);
            return (offspring);
        }//cross1 ends

        public byte[][] cross_over2(byte[][] inda, byte[][] indb, float just) {
            int row = inda.length; //# rows
            int col = inda[0].length; //# columns in row 0
            int rand;
            byte[][] offspring = new byte[row][col];
            for (int i = 0; i < row; i++) {
                for (int j = 0; j < col; j++) {
                    offspring[i][j] = inda[i][j];
                }
            }
//choose a random row from both parents
            int randa = rand_int(row);
            int randb = rand_int(row);
//recipe of row randa in offspring, if value of template[i]=1
//then offspring[randa][i]=indb[randb][i], if template[i]=0
//then offspring[randa][i]=inda[randa][i]
            int[] template = new int[col];
            for (int i = 0; i < col; i++) {
                rand = rand_int(100);
                if (rand > just - 1) {
                    template[i] = 1;
                }
            }//for ends
            for (int i = 0; i < col; i++) {
                if (template[i] == 1) {
                    offspring[randa][i] = indb[randb][i];
                }
            }
            offspring = correct_Matrix(offspring, randa);
            return (offspring);
        }//cross2 ends

        public void printMatrix(byte[][] matrix) {
            int row = matrix.length, col = matrix[0].length;
            System.out.println();
            for (int i = 0; i < row; i++) {
                for (int j = 0; j < col; j++) {
                    System.out.print(matrix[i][j]);
                }
                System.out.println();
            }
        }
//random customer chosen and moved to another car

        public void mutation(byte[][] matrix) {
            int i = 0, out = 0, row = matrix.length, col = matrix[0].length;
            int randcust = rand_int(col);
            while (out == 0) {
                out = matrix[i][randcust];
                i++;
            }
            matrix[i - 1][randcust] = 0;
            int rand = rand_int(row);
            while (rand == i - 1) {
                rand = rand_int(row);
            }
            matrix[rand][randcust] = 1;
        }
//returns a random integer on the interval [0, L-1]

        public int rand_int(int L) {
            return ((int) (Math.random() * L));
        }

        public byte[][] correct_Matrix(byte[][] matrix, int randa) {
            int rand, row = matrix.length, col = matrix[0].length;
            for (int j = 0; j < row; j++) {
                matrix[j][0] = 1;
            }
            for (int i = 1; i < col; i++) {
                int su = 0;
                for (int j = 0; j < row; j++) {
                    su = su + matrix[j][i];
                }
                switch (su) {
                    case 0:
                        rand = rand_int(row);
                        while (rand == randa) {
                            rand = rand_int(row);
                        }
                        matrix[rand][i] = 1;
                        break;
                    case 1:
                        break;
                    case 2:
                        for (int j = 0; j < row; j++) {
                            if (matrix[j][i] == 1 && j != randa) {
                                matrix[j][i] = 0;
                            }
                        }
                        break;
                    default:
                        System.out.println("Something has gone wrong in cross_over");
                }//switch ends
            }//for ends
            return (matrix);
        }//correct_Matrix ends

        public int getUTimeWindow(int c) {
            return (super.getUTimeWindow(c));
        }
    }//class ends

}
