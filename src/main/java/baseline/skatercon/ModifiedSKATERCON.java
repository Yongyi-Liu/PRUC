package baseline.skatercon;

import util.Area;
import util.Region;
import baseline.skater.ModifiedSKATER;

import java.util.ArrayList;

/**
 * This class implements the Modified SKATERCON as one of our baseline algorithms
 * The original SKATER is presented in Aydin, Orhun, Mark V. Janikas, Renato Assunçao, and Ting-Hwan Lee. "SKATER-CON: Unsupervised regionalization via stochastic tree partitioning within a consensus framework using random spanning trees." In Proceedings of the 2nd ACM SIGSPATIAL International Workshop on AI for Geographic Knowledge Discovery, pp. 33-42. 2018.
 */
public class ModifiedSKATERCON {
    private int thread_num;
    Region[][] skater_results;
    Region[] final_regions;
    private long runtime;
    private ArrayList<Area> all_areas;
    private int coarsen_threshold;
    private int sc;
    private int p;
    private long threshold;

    /**
     *
     * @param all_areas The input areas
     * @param sc The sc parameter in SKATER
     * @param p The number of regions
     * @param threshold The threshold on the user-defined value, when set to 0, the Modified SAKTERCON becomes SKATERCON
     */
    public ModifiedSKATERCON(ArrayList<Area> all_areas, int sc , int p , long threshold) throws InterruptedException, CloneNotSupportedException {
        long start_runtime = System.currentTimeMillis();
        this.thread_num = 4; //the number of random spanning tree as input
        this.all_areas = all_areas;
        this.skater_results = new Region[thread_num][];
        this.coarsen_threshold = 100;
        this.sc = sc;
        this.p = p;
        this.threshold = threshold;
        skater_con_run();
        long end_runtime = System.currentTimeMillis();
        runtime = end_runtime - start_runtime;
    }

    public long getRuntime()
    {
        return runtime;
    }


    private void skater_con_run() throws CloneNotSupportedException, InterruptedException {
        Parallel_SKATERCON[] skater_threads = new Parallel_SKATERCON[thread_num];


        for(int i = 0 ; i < thread_num ; i++)
        {
            skater_threads[i] = new Parallel_SKATERCON(i , all_areas , sc , p , skater_results);
            skater_threads[i].start();
        }

        for(int i = 0 ; i < thread_num ; i++)
        {
            skater_threads[i].join();
        }


        ArrayList<Region[]> regions_sets = new ArrayList<>();

        for(Region[] regions : skater_results)
        {
            if(regions != null)
            {
                regions_sets.add(regions);
            }
        }



        if(regions_sets.size() == 1)
        {
            final_regions = regions_sets.get(0);
            return;
        }

        if(regions_sets.size() == 0)
        {
            final_regions = null;
            return;
        }


        this.final_regions = new Metis(all_areas , regions_sets , p , coarsen_threshold , threshold).Metis_start();
    }

    public Region[] getFinal_regions() {
        return final_regions;
    }

    class Parallel_SKATERCON extends Thread
    {
        ArrayList<Area> all_areas_local;
        int sc;
        int index;
        int p;
        Region[][] all_regions_arr;

        public Parallel_SKATERCON(int index , ArrayList<Area> all_areas , int sc , int p , Region[][] all_regions_arr) throws CloneNotSupportedException {
            this.index = index;
            this.sc = sc;
            this.p = p;
            this.all_areas_local = Area.area_list_copy(all_areas);
            this.all_regions_arr = all_regions_arr;
        }

        public void run()
        {

            Region[] result_regions;
            ArrayList<Area> rst = new RandomSpanningTree(all_areas_local).run_RST();

            result_regions = new ModifiedSKATER(rst,  sc , p ,  threshold, true).get_regions();
            all_regions_arr[index] = result_regions;

        }
    }









}




