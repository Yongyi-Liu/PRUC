package test;


import GSLO.LocalOptimization;
import GSLO.GlobalSearch;
import baseline.greedy.*;
import baseline.skater.ModifiedSKATER;
import baseline.skatercon.ModifiedSKATERCON;
import util.Area;
import util.Preprocess;
import util.Region;

import java.io.IOException;
import java.util.ArrayList;

public class TestGeneral {


    public static void main(String args[]) throws InterruptedException, CloneNotSupportedException, IOException {

        int GSLO_iter = 100;
        int iter_skatercon = 100;
        int iter_prob = 100;
        double threshold_scale = 0.02;
        int default_p = 10;


        var_iter_seed_quality(GSLO_iter , threshold_scale , default_p  , "2k");

        var_max_no_move(GSLO_iter , threshold_scale ,  "2k");

        test_time_breakdown(GSLO_iter , default_p , threshold_scale , "2k" );

        greedy_baseline(GSLO_iter  , threshold_scale , "2k" , false);

        greedy_baseline(GSLO_iter , threshold_scale , "island" , true);

        var_p_measure_hetero_runtime(GSLO_iter , true , iter_skatercon , threshold_scale , "2k");
        var_p_measure_hetero_runtime(GSLO_iter , true , iter_skatercon , threshold_scale , "diversity");

        var_thre_measure_hetero_runtime(GSLO_iter , true , iter_skatercon , default_p , threshold_scale , "2k");
        var_thre_measure_hetero_runtime(GSLO_iter , true , iter_skatercon , default_p , threshold_scale , "diversity");


        var_p_measure_hetero_runtime(GSLO_iter , true , iter_skatercon , 0 , "2k");
        var_p_measure_hetero_runtime(GSLO_iter , true , iter_skatercon , 0 , "diversity");


        var_thre_finding_feasible_partition(iter_prob , threshold_scale , "2k");
        var_thre_finding_feasible_partition(iter_prob , threshold_scale , "diversity");
        var_p_finding_feasible_partition(iter_prob , threshold_scale , "2k");
        var_p_finding_feasible_partition(iter_prob , threshold_scale , "diversity");

        var_dataset_hetero_runtime(GSLO_iter , default_p , true , iter_skatercon,  threshold_scale);

    }


    /**
     * This method tests the time and heterogeneity of GSLO, SKATER, and SKATERCON among different datasets
     * This corresponds to Experiment in Section 7.2.8
     * @param GSLO_iter the number of iterations for GSLO
     * @param p the number of predefined regions
     * @param run_skater whether or not to run SKATER
     * @param iter_skatercon the number of iterations in SKATERCON
     * @param scale
     */
    public static void var_dataset_hetero_runtime(int GSLO_iter , int p , boolean run_skater , int iter_skatercon , double scale) throws IOException, CloneNotSupportedException, InterruptedException {

        System.out.println(" Experiment starts, this experiment tests the runtime and heterogeneity under different size datasets");
        GlobalSearch sol;
        String[] datasets = new String[]{"2k" , "5k" , "10k" ,"20k" , "30k","40k"};
        for(String dataset : datasets)
        {
            System.out.println("the current dataset is the " + dataset);
            ArrayList<Area> all_areas = Preprocess.GeoSetBuilder(dataset);
            long total_ext = 0;
            for(Area area : all_areas)
            {
                total_ext += area.get_extensive_attr();
            }
            long threshold = (long)(scale * total_ext);
            System.out.println("the total extensive attribtue is " + total_ext + " the scaled threshold is " + threshold);

            ArrayList<Long> GSLO_hetero = new ArrayList<>();
            ArrayList<Long> GSLO_runtime = new ArrayList<>();
            ArrayList<Long> SKATERCON_hetero = new ArrayList<>();
            ArrayList<Long> SKATERCON_runtime = new ArrayList<>();

            for(int i = 0 ; i < GSLO_iter ; i++)
            {

                sol = new GlobalSearch(Area.area_list_copy(all_areas) , p , all_areas.size() , threshold , false);
                if(sol.solved())
                {
                    LocalOptimization lo = new LocalOptimization(sol , all_areas.size() , 0.99 , sol.get_all_areas() , sol.get_regions() , threshold);
                    Region.test_result_correctness(sol.get_regions() , all_areas , threshold , true);
                    GSLO_hetero.add(lo.getBest_hetero());
                    GSLO_runtime.add(sol.getTotal_running_time() + lo.getTotal_time());
                }

                else {
                    GSLO_runtime.add(sol.getTotal_running_time());
                }

            }

            System.out.println("GSLO hetero" + compute_long_ave(GSLO_hetero));
            System.out.println("GSLO runtime " + compute_long_ave(GSLO_runtime) );

            if(run_skater)
            {

                ModifiedSKATER m_skater = new ModifiedSKATER(Area.area_list_copy(all_areas) , 30 ,  p , threshold, false);
                if(m_skater.getTrees() == null)
                {
                    System.out.println("skater did not find a solution with p = " + p + " the runtime is " + m_skater.getRuntime());
                }

                else
                {
                    Region.test_result_correctness(m_skater.get_regions() , all_areas , threshold , false);
                    System.out.println("skater hetero is "+ Region.get_all_region_hetero(m_skater.get_regions()) + " runtime " + m_skater.getRuntime());
                }
            }


            for(int i = 0; i < iter_skatercon; i++)
            {
                ModifiedSKATERCON mskc = new ModifiedSKATERCON(Area.area_list_copy(all_areas) , 30 , p , threshold);
                SKATERCON_runtime.add(mskc.getRuntime());
                if(mskc.getFinal_regions() != null)
                {
                    Region.test_result_correctness(mskc.getFinal_regions() , all_areas , threshold , false);
                    SKATERCON_hetero.add(Region.get_all_region_hetero(mskc.getFinal_regions()));
                }
            }
            System.out.println("skatercon hetero is " + compute_long_ave(SKATERCON_hetero) + " runtime " + compute_long_ave(SKATERCON_runtime));
            System.out.println();
        }

    }


    /**
     * This method tests the heterogeneity and runtime under different p under the default dataset.
     * This corresponds to the experiment in Section 7.2.4
     * If threshold = 0, then this corresponds to the p-regions test in Section 7.2.7
     * @param iter_GSLO the number of iterations of GSLO
     * @param run_skater whether or not to run SKATER
     * @param iter_time_skatercon the number of iterations in SKATERCON
     * @param thre_scale the threshold scale, i.e., the threshold on each region = total extensive attribute * thre_scale
     * @param dataset the dataset to compute
     */
    public static void var_p_measure_hetero_runtime(int iter_GSLO , boolean run_skater , int iter_time_skatercon , double thre_scale , String dataset) throws IOException, CloneNotSupportedException, InterruptedException {
        ArrayList<Area> all_areas = Preprocess.GeoSetBuilder(dataset);
        long base_ext = 0;
        for(Area area : all_areas)
        {
            base_ext += area.get_extensive_attr();
        }
        base_ext = (long)(thre_scale * base_ext);

        System.out.println("var_p_measure_hetero_runtime,  the performance measure is heterogeneity and runtime" + " the dataset is " + dataset);
        if(base_ext == 0)
        {
            System.out.println("this is p-regions test");
        }

        long threshold = base_ext;
        GlobalSearch sol;
        for(int p = 5 ; p <= 50 ; p = p + 5)
        {
            System.out.println("the current p is" + p);
            ArrayList<Long> GSLO_hetero = new ArrayList<>();
            ArrayList<Long> GSLO_runtime = new ArrayList<>();
            ArrayList<Long> skater_con_results = new ArrayList<>();
            ArrayList<Long> skater_con_runtime = new ArrayList<>();

            for(int i = 0 ; i < iter_GSLO ; i++)
            {

                sol = new GlobalSearch(Area.area_list_copy(all_areas) , p , all_areas.size(), threshold , false);
                if(sol.solved())
                {
                    LocalOptimization lo = new LocalOptimization(sol , all_areas.size() , 0.99 , sol.get_all_areas() , sol.get_regions() , threshold);
                    GSLO_hetero.add(lo.getBest_hetero());
                    GSLO_runtime.add(sol.getTotal_running_time() + lo.getTotal_time());
                    Region.test_result_correctness(sol.get_regions() , all_areas , threshold, true);
                }

                else
                {
                    GSLO_runtime.add(sol.getTotal_running_time());
                }



            }

            System.out.println(" GSLO hetero " + compute_long_ave(GSLO_hetero));
            System.out.println(" GSLO runtime " + compute_long_ave(GSLO_runtime));

            if(run_skater)
            {

                ModifiedSKATER m_skater = new ModifiedSKATER(Area.area_list_copy(all_areas) , 30 ,  p , threshold , false);
                if(m_skater.getTrees() == null)
                {
                    System.out.println("skater did not find a solution with p = " + p + " the runtime is " + m_skater.getRuntime());
                }

                else
                {
                    Region.test_result_correctness(m_skater.get_regions() , all_areas , threshold , false);
                    System.out.println("skater hetero is " + Region.get_all_region_hetero(m_skater.get_regions()) + " runtime " + m_skater.getRuntime());
                }
            }


            for(int i = 0; i < iter_time_skatercon; i++)
            {
                ModifiedSKATERCON mskc = new ModifiedSKATERCON(Area.area_list_copy(all_areas) , 30 , p , threshold);
                skater_con_runtime.add(mskc.getRuntime());
                if(mskc.getFinal_regions() != null)
                {
                    Region.test_result_correctness(mskc.getFinal_regions() , all_areas , threshold , false);
                    skater_con_results.add(Region.get_all_region_hetero(mskc.getFinal_regions()));
                }
            }
            System.out.println("heterogeneity skatercon " + compute_long_ave(skater_con_results) + " runtime " + compute_long_ave(skater_con_runtime));
            System.out.println();
        }

    }

    /**
     * This method tests the heterogeneity and runtime under different threshold under the default dataset.
     * This corresponds to the experiment in Section 7.2.5
     * @param iter_GSLO the number of iterations in GSLO
     * @param run_skater whether or not to run SKATER
     * @param iter_time_skatercon the number of iterations in SKATERCON
     * @param p the number of predefined regions
     * @param thre_scale the threshold scale, i.e., the threshold on each region = total extensive attribute * thre_scale
     * @param dataset the dataset to compute

     */
    public static void var_thre_measure_hetero_runtime(int iter_GSLO , boolean run_skater , int iter_time_skatercon , int p , double thre_scale,  String dataset) throws IOException, CloneNotSupportedException, InterruptedException {
        ArrayList<Area> all_areas = Preprocess.GeoSetBuilder(dataset);
        System.out.println("Experiment starts, the variable is threshold, the performance measure is heterogeneity and runtime");
        System.out.println("the current dataset is the " + dataset);
        GlobalSearch sol;


        long all_ext = 0;
        for(Area area : all_areas)
        {
            all_ext += area.get_extensive_attr();
        }
        long starting_threshold = (long)(all_ext * 0.5 * thre_scale);
        System.out.println("the total extensive attribute is " + all_ext + " the starting threshold is " + starting_threshold);

        for(long threshold = starting_threshold; threshold <= 10 * starting_threshold ; threshold += starting_threshold)
        {
            System.out.println("the current threshold is" + threshold);
            ArrayList<Long> GSLO_hetero = new ArrayList<>();
            ArrayList<Long> GSLO_runtime = new ArrayList<>();
            ArrayList<Long> skatercon_hetero = new ArrayList<>();
            ArrayList<Long> skatercon_runtime = new ArrayList<>();

            for(int i = 0 ; i < iter_GSLO ; i++)
            {

                sol = new GlobalSearch(Area.area_list_copy(all_areas) , p , all_areas.size(), threshold , false);
                if(sol.solved())
                {
                    LocalOptimization up = new LocalOptimization(sol , all_areas.size() , 0.99 , sol.get_all_areas() , sol.get_regions() , threshold);
                    GSLO_hetero.add(up.getBest_hetero());
                    GSLO_runtime.add(sol.getTotal_running_time() + up.getTotal_time());

                    Region.test_result_correctness(sol.get_regions() , all_areas , threshold, true);
                }

                else
                {
                    GSLO_runtime.add(sol.getTotal_running_time());
                }



            }


            System.out.println("GSLO hetero " + compute_long_ave(GSLO_hetero));
            System.out.println("GSLO runtime " + compute_long_ave(GSLO_runtime));


            if(run_skater)
            {
                ModifiedSKATER m_skater = new ModifiedSKATER(Area.area_list_copy(all_areas) , 30 ,  p , threshold, false);


                if(m_skater.getTrees() == null)
                {
                    System.out.println("skater did not find a solution with p = " + p + " the runtime is " + m_skater.getRuntime());
                }

                else
                {
                    Region.test_result_correctness(m_skater.get_regions() , all_areas , threshold , false);
                    System.out.println("skater hetero " + Region.get_all_region_hetero(m_skater.get_regions()) + " runtime " + m_skater.getRuntime());
                }
            }

            for(int i = 0; i < iter_time_skatercon; i++)
            {
                ModifiedSKATERCON mskc = new ModifiedSKATERCON(Area.area_list_copy(all_areas) , 30 , p , threshold);
                skatercon_runtime.add(mskc.getRuntime());
                if(mskc.getFinal_regions() != null)
                {
                    Region.test_result_correctness(mskc.getFinal_regions() , all_areas , threshold , false);
                    skatercon_hetero.add(Region.get_all_region_hetero(mskc.getFinal_regions()));
                }
            }

            System.out.println("heterogeneity skatercon " + compute_long_ave(skatercon_hetero) + " runtime " + compute_long_ave(skatercon_runtime));
            System.out.println();

        }

    }

    /**
     * This method computes the effectiveness, i.e, the probability of finding a feasible partition, of GSLO, SKATER, and SKATERCON under different p
     * This corresponds to the Experiment in Section 7.2.6
     * @param iter_num the number of iterations for GSLO and skatercon
     * @param thre_scale the threshold scale, i.e., the threshold on each region = total extensive attribute * thre_scale
     * @param dataset the dataset to compute
     */
    public static void var_p_finding_feasible_partition(int iter_num, double thre_scale, String dataset) throws IOException, CloneNotSupportedException, InterruptedException {
        ArrayList<Area> all_areas = Preprocess.GeoSetBuilder(dataset);
        System.out.println("var_p_finding_feasible_partition starts " + " the dataset is " + dataset);

        long total_ext = 0;

        for(Area area : all_areas)
        {
            total_ext += area.get_extensive_attr();
        }

        long threshold = (long)(total_ext * thre_scale);

        for(int p = 5 ; p <= 50 ; p = p+5)
        {
            int GSLO_fail = 0;
            int skater_con_fail = 0;
            GlobalSearch sol;
            for(int i = 0 ; i < iter_num ; i++)
            {

                sol = new GlobalSearch(Area.area_list_copy(all_areas) , p , all_areas.size(), threshold , false);
                if(!sol.solved()) {
                    GSLO_fail += 1;
                }


                ModifiedSKATERCON mskc = new ModifiedSKATERCON(Area.area_list_copy(all_areas) , 30 , p , threshold);

                if(mskc.getFinal_regions() == null)
                {
                    skater_con_fail += 1;
                }
            }

            ModifiedSKATER m_skater = new ModifiedSKATER(Area.area_list_copy(all_areas) , 30 ,  p , threshold, false);

            System.out.println("skater failed " + (m_skater.getTrees() == null));
            System.out.println("the current p is "+ p + " the current threshold is " + threshold +  " GSLO fail " + GSLO_fail + "skatercon failure" +  skater_con_fail);

        }
    }


    /**
     * This method computes the effectiveness, i.e, the probability of finding a feasible partition, of GSLO, SKATER, and SKATERCON under different threshold
     * This corresponds to the Experiment in Section 7.2.6
     * @param iter_num the number of iterations for GSLO and skatercon
     * @param thre_scale the threshold scale, i.e., the threshold on each region = total extensive attribute * thre_scale
     * @param dataset the dataset to compute
     */
    public static void var_thre_finding_feasible_partition(int iter_num, double thre_scale, String dataset) throws IOException, CloneNotSupportedException, InterruptedException {
        ArrayList<Area> all_areas = Preprocess.GeoSetBuilder(dataset);
        System.out.println("var_thre_finding_feasible_partition starts " + " the dataset is " + dataset);
        long total_ex = 0;
        for(Area area : all_areas)
        {
            total_ex += area.get_extensive_attr();
        }
        long threshold = (long)(0.5 * total_ex * thre_scale);
        long base_thre = threshold;
        long aim_threshold = threshold * 10;
        int p = 10;
        GlobalSearch sol;
        for(; threshold < aim_threshold ; threshold += base_thre)
        {
            int pruc_hg_fail = 0;
            int pruc_rg_fail = 0;
            int skater_con_fail = 0;

            for(int i = 0 ; i < iter_num ; i++)
            {

                sol = new GlobalSearch(Area.area_list_copy(all_areas) , p , all_areas.size(), threshold , false);
                if(!sol.solved()) {
                    pruc_rg_fail += 1;
                }

                ModifiedSKATERCON mskc = new ModifiedSKATERCON(Area.area_list_copy(all_areas) , 30 , p , threshold);

                if(mskc.getFinal_regions() == null)
                {
                    skater_con_fail += 1;
                }
            }

            ModifiedSKATER m_skater = new ModifiedSKATER(Area.area_list_copy(all_areas) , 30 ,  p , threshold, false);


            System.out.println("skater failed?" + (m_skater.getTrees() == null));
            System.out.println("the current p is "+ p + " the current threshold is " + threshold + " pruc_hg failure " + pruc_hg_fail + " pruc_rg_failure" + pruc_rg_fail + "skatercon failure" +  skater_con_fail);

        }
    }

    /**
     * This method tests how the number of iterations in Local Optimization affect the heterogeneity and runtime of GSLO
     * This corresponds to the experiment in Section 7.1 GSLO Parameter Tuning Part 2
     * @param GSLO_iter the number of iterations for GSLO
     * @param thre_scale the threshold scale, i.e., the threshold on each region = total extensive attribute * thre_scale
     * @param dataset the dataset to compute
     */
    public static void var_max_no_move(int GSLO_iter, double thre_scale, String dataset) throws IOException, CloneNotSupportedException, InterruptedException {
        System.out.println("var_max_no_move starts" + " the dataset is " + dataset);
        ArrayList<Area> all_areas = Preprocess.GeoSetBuilder(dataset);

        long total_ext = 0;
        for(Area area : all_areas)
        {
            total_ext += area.get_extensive_attr();
        }
        System.out.println("the total extensive attribute is " + total_ext);

        int p = 10;
        ArrayList<Long> reduction0001 = new ArrayList<>();
        ArrayList<Long> time0001 = new ArrayList<>();

        ArrayList<Long> reduction001 = new ArrayList<>();
        ArrayList<Long> time001 = new ArrayList<>();

        ArrayList<Long> reduction01 = new ArrayList<>();
        ArrayList<Long> time01 = new ArrayList<>();

        ArrayList<Long> reduction1 = new ArrayList<>();
        ArrayList<Long> time1 = new ArrayList<>();

        ArrayList<Long> reduction10 = new ArrayList<>();
        ArrayList<Long> time10 = new ArrayList<>();

        ArrayList<Long> reduction100 = new ArrayList<>();
        ArrayList<Long> time100 = new ArrayList<>();






        int case0001 = (int)(0.001 * all_areas.size());
        int case001 = (int)(0.01 * all_areas.size());
        int case01 = (int)(0.1 * all_areas.size());
        int case1 = 1 * all_areas.size();
        int case10 = 10 * all_areas.size();
        int case100 = 100 * all_areas.size();



        long threshold = (long)(total_ext * thre_scale);
        System.out.println("the scaled total extensive attribute is " + threshold);
        GlobalSearch sol;

        for(int i = 0 ; i < GSLO_iter ; i ++)
        {
            sol = new GlobalSearch(Area.area_list_copy(all_areas) , p , all_areas.size(), threshold , false);
            LocalOptimization up;
            if(sol.solved())
            {
                Region.test_result_correctness(sol.get_regions() , all_areas , threshold , true);
                up = new LocalOptimization(sol , case0001, 0.99 , sol.get_all_areas() , sol.get_regions() , threshold);
                time0001.add(up.getTotal_time() + sol.getTotal_running_time());
                reduction0001.add(up.getBest_hetero());
            }


            sol = new GlobalSearch(Area.area_list_copy(all_areas) , p , all_areas.size(), threshold , false);
            if(sol.solved())
            {
                up = new LocalOptimization(sol , case001, 0.99 , sol.get_all_areas() , sol.get_regions()  , threshold);
                time001.add(up.getTotal_time() + sol.getTotal_running_time());
                reduction001.add(up.getBest_hetero());
            }

            sol = new GlobalSearch(Area.area_list_copy(all_areas) , p , all_areas.size(), threshold , false);
            if(sol.solved())
            {
                up = new LocalOptimization(sol , case01, 0.99 ,sol.get_all_areas() , sol.get_regions() , threshold);
                time01.add(up.getTotal_time() + sol.getTotal_running_time());
                reduction01.add(up.getBest_hetero());
            }

            sol = new GlobalSearch(Area.area_list_copy(all_areas) , p , all_areas.size(), threshold , false);
            if(sol.solved())
            {
                up = new LocalOptimization(sol , case1, 0.99 , sol.get_all_areas() , sol.get_regions() , threshold);
                time1.add(up.getTotal_time() + sol.getTotal_running_time());
                reduction1.add(up.getBest_hetero());
            }

            sol = new GlobalSearch(Area.area_list_copy(all_areas) , p , all_areas.size(), threshold, false);
            if(sol.solved())
            {
                up = new LocalOptimization(sol , case10, 0.99 , sol.get_all_areas() , sol.get_regions() , threshold);
                time10.add(up.getTotal_time() + sol.getTotal_running_time());
                reduction10.add(up.getBest_hetero());
            }

            sol = new GlobalSearch(Area.area_list_copy(all_areas) , p , all_areas.size(), threshold , false);
            if(sol.solved())
            {
                up = new LocalOptimization(sol , case100, 0.99 , sol.get_all_areas() , sol.get_regions() , threshold);
                time100.add(up.getTotal_time() + sol.getTotal_running_time());
                reduction100.add(up.getBest_hetero());
            }

        }


        System.out.println("when max_no_improve = 0.001 * size " + "the avg improve is " + compute_long_ave(reduction0001) + " avg time is " + compute_long_ave(time0001));
        System.out.println("when max_no_improve = 0.01 * size " + "the avg improve is " + compute_long_ave(reduction001) + " avg time is " + compute_long_ave(time001));
        System.out.println("when max_no_improve = 0.1 * size " + "the avg improve is " + compute_long_ave(reduction01) + " avg time is " + compute_long_ave(time01));
        System.out.println("when max_no_improve = 1 * size " + "the avg improve is " + compute_long_ave(reduction1) + " avg time is " + compute_long_ave(time1));
        System.out.println("when max_no_improve = 10 * size " + "the avg improve is " + compute_long_ave(reduction10) + " avg time is " + compute_long_ave(time10));
        System.out.println("when max_no_improve = 100 * size " + "the avg improve is " + compute_long_ave(reduction100) + " avg time is " + compute_long_ave(time100));

    }

    /**
     * The method tests how the number of iterations in Seed Identification affect the seed quality, heterogeneity, and runtime
     * This corresponds to Section 7.1 GSLO Parameter Tuning Part 1
     * @param GSLO_iter the number of iterations of GSLO
     * @param thre_scale the threshold scale, i.e., the threshold on each region = total extensive attribute * thre_scale
     * @param p the predefined number of regions
     * @param dataset the dataset to compute
     */
    public static void var_iter_seed_quality(int GSLO_iter, double thre_scale , int p , String dataset) throws IOException, CloneNotSupportedException, InterruptedException {
        System.out.println(" test how the number iterations in seeding affect the  " + " the dataset is " + dataset);
        ArrayList<Area> all_areas = Preprocess.GeoSetBuilder(dataset);
        long total_ex = 0;
        for(Area area : all_areas)
        {
            total_ex += area.get_extensive_attr();
        }
        System.out.println();

        long threshold = (long)(total_ex * thre_scale);
        System.out.println("the total extensive attribute is " + total_ex + " the computed threshold is " + threshold);

        int datasetsize = all_areas.size();

        //int[] parameters = new int[]{};
        int[] parameters = new int[]{0 , (int)(0.001 * datasetsize) , (int)(0.01 * datasetsize) , (int)(0.1 * datasetsize) , datasetsize , 10 * datasetsize , 100 * datasetsize , 1000 * datasetsize};

        GlobalSearch sol;
        for(int parameter : parameters)
        {
            ArrayList<Long> hetero = new ArrayList<>();
            ArrayList<Long> time = new ArrayList<>();
            ArrayList<Double> seed_quality = new ArrayList<>();
            for(int i = 0 ; i < GSLO_iter ; i++)
            {
                sol = new GlobalSearch(Area.area_list_copy(all_areas) , p , parameter, threshold , false);
                LocalOptimization uh = new LocalOptimization(sol , all_areas.size() , 0.99 , sol.get_all_areas() , sol.get_regions() , threshold);
                if(sol.solved())
                {
                    Region.test_result_correctness(sol.get_regions() , all_areas , threshold , true);
                    hetero.add(uh.getBest_hetero());
                    time.add(sol.getTotal_running_time() + uh.getTotal_time());
                    seed_quality.add(sol.get_seed_quality());
                }

            }
            System.out.println(" the current scale size is "  + parameter + " avg hetero " + compute_long_ave(hetero) + " avg seed quality " + compute_double_ave(seed_quality)  + " avg runtime " + compute_long_ave(time));

        }


        System.out.println(" the following test the k-means++ seeding ");
        ArrayList<Long> hetero = new ArrayList<>();
        ArrayList<Long> time = new ArrayList<>();
        for(int i = 0 ; i < GSLO_iter ; i++)
        {
            sol = new GlobalSearch(Area.area_list_copy(all_areas) , p , -1, threshold , false);
            LocalOptimization uh = new LocalOptimization(sol , all_areas.size() , 0.99 , sol.get_all_areas() , sol.get_regions() , threshold);
            if(sol.solved())
            {
                Region.test_result_correctness(sol.get_regions() , all_areas , threshold , true);
                hetero.add(uh.getBest_hetero());
                time.add(sol.getTotal_running_time() + uh.getTotal_time());
            }
        }
        System.out.println(" k-means++ seeding "  + " avg hetero " + compute_long_ave(hetero)  + " avg runtime " + compute_long_ave(time));




    }


    /**
     * This method compares GSLO with the greedy baseline
     * This corresponds to Section 7.2.2 and 7.2.3
     * @param iternum The number of iterations for GSLO and simple greedy
     * @param thre_scale the threshold scale, i.e., the threshold on each region = total extensive attribute * thre_scale
     * @param dataset the datset to compute
     * @param island whether or not to consider island
     * @throws IOException
     * @throws CloneNotSupportedException
     * @throws InterruptedException
     */
    public static void greedy_baseline(int iternum ,  double thre_scale , String dataset , boolean island) throws IOException, CloneNotSupportedException, InterruptedException
    {
        System.out.println("compare GSLO against simple greedy " + " the dataset is " + dataset + " island area present? " + island);
        ArrayList<Area> all_areas = Preprocess.GeoSetBuilder(dataset);
        long total_ex = 0;
        for(Area area : all_areas)
        {
            total_ex += area.get_extensive_attr();
        }
        System.out.println();

        long threshold = (long)(total_ex * thre_scale);

        for(int p = 5; p <= 50 ; p += 5)
        {
            System.out.println("the current p is " + p);
            ArrayList<Long> pruc_hetero = new ArrayList<>();
            ArrayList<Long> pruc_runtime = new ArrayList<>();
            int pruc_success = 0;

            ArrayList<Long> sg_hetero = new ArrayList<>();
            ArrayList<Long> sg_runtime = new ArrayList<>();
            int sg_success = 0;

            for(int i = 0 ; i < iternum ; i++)
            {
                GlobalSearch sol = new GlobalSearch(Area.area_list_copy(all_areas) , p, all_areas.size(), threshold , island);
                if(sol.solved())
                {
                    LocalOptimization up = new LocalOptimization(sol , all_areas.size() , 0.99 , sol.get_all_areas() , sol.get_regions() , threshold);
                    pruc_hetero.add(up.getBest_hetero());
                    pruc_runtime.add(sol.getTotal_running_time() + up.getTotal_time());
                    pruc_success += 1;
                }

                else
                {
                    pruc_runtime.add(sol.getTotal_running_time());
                }

                SimpleGreedy sg = new SimpleGreedy(Area.area_list_copy(all_areas) , p , threshold , island);
                if(sg.get_solved())
                {
                    GreedySearch hb = new GreedySearch(sg , sg.getAll_areas() , sg.getRegions() , threshold);
                    sg_hetero.add(hb.getBest_hetero());
                    sg_runtime.add(hb.getTotal_time() + sg.getTotal_running_time());
                    sg_success += 1;

                }
                else
                {
                    sg_runtime.add(sg.getTotal_running_time());
                }


            }

            System.out.println("pruc heterogeneity " + compute_long_ave(pruc_hetero) + " runtime " + compute_long_ave(pruc_runtime) + " success rate " + pruc_success);
            System.out.println("sg   heterogeneity " + compute_long_ave(sg_hetero) + " runtime " + compute_long_ave(sg_runtime) + " success rate " + sg_success );

        }
    }

    /**
     * This method tests the runtime on different phases in GSLO
     * This corresponds to the experiment in Section 7.2.1 Time Breakdown Analysis
     * @param GSLO_iter the number of iterations in GSLO
     * @param default_p the default value on the predefined number of regions
     * @param thre_scale the threshold scale, i.e., the threshold on each region = total extensive attribute * thre_scale
     * @param dataset the dataset to compute
     */
    public static void test_time_breakdown(int GSLO_iter, int default_p, double thre_scale , String dataset) throws IOException, CloneNotSupportedException, InterruptedException
    {
        System.out.println(" testing time breakdown, the time of different phases ");
        ArrayList<Area> all_areas = Preprocess.GeoSetBuilder(dataset);
        long all_ext = 0;
        for(Area area : all_areas)
        {
            all_ext += area.get_extensive_attr();
        }
        System.out.println();

        long starting_threshold = (long)(all_ext * 0.5 * thre_scale);
        System.out.println("changing the threshold");
        for(long threshold = starting_threshold * 8; threshold <= 10 * starting_threshold ; threshold += starting_threshold)
        {

            ArrayList<Long> seed_times = new ArrayList<>();
            ArrayList<Long> grow_times = new ArrayList<>();
            ArrayList<Long> enclaves_times = new ArrayList<>();
            ArrayList<Long> interregion_times = new ArrayList<>();
            ArrayList<Long> flow_times = new ArrayList<>();
            ArrayList<Long> before_heteros = new ArrayList<>();
            ArrayList<Long> local_times = new ArrayList<>();
            ArrayList<Long> after_heteros = new ArrayList<>();
            ArrayList<Double> rates = new ArrayList<>();
            int solved_cases = 0;

            for(int i = 0 ; i < GSLO_iter ; i++)
            {
                GlobalSearch sol = new GlobalSearch(Area.area_list_copy(all_areas) , default_p , all_areas.size(), threshold , false);
                seed_times.add(sol.getSeed_time());
                grow_times.add(sol.getRegion_growth_time());
                enclaves_times.add(sol.getEnclaves_assign_time());

                if(sol.isInterregion_flag())
                {
                    interregion_times.add(sol.getInterregion_update_time());
                }

                if(sol.isFlow_flag())
                {
                    flow_times.add(sol.getIndirect_flow_time());
                }

                if(sol.solved())
                {
                    solved_cases ++;
                    long before_hetero = Region.get_all_region_hetero(sol.get_regions());
                    LocalOptimization up = new LocalOptimization(sol , all_areas.size() , 0.99 , sol.get_all_areas() , sol.get_regions() , threshold);
                    local_times.add(up.getTotal_time());
                    long after_hetero = up.getBest_hetero();
                    double improve = ((before_hetero - after_hetero) * 1.0) / before_hetero;
                    before_heteros.add(before_hetero);
                    after_heteros.add(after_hetero);
                    rates.add(improve);

                }




            }

            System.out.println("the current threshold is " + threshold + " the seed runtime is " + compute_long_ave(seed_times) + " region growth runtime " + compute_long_ave(grow_times) + " enclaves time " + compute_long_ave(enclaves_times));
            if(interregion_times.size() > 0)
            {
                System.out.println("num time entered interregion " + interregion_times.size() + "the interregional time is " + compute_long_ave(interregion_times));
            }

            if(flow_times.size() > 0)
            {
                System.out.println("num times entered flow " + flow_times.size() +   " the flow time is " + compute_long_ave(flow_times));
            }

            if(after_heteros.size() > 0)
            {
                System.out.println("avg before is " + compute_long_ave(before_heteros) + " ave after is " + compute_long_ave(after_heteros) + "rate is " + compute_double_ave(rates) + " local time is " + compute_long_ave(local_times));
            }

            System.out.println("the total number of solved cases is " + solved_cases);

        }

        System.out.println("var p");
        long threshold = (long)(all_ext  * thre_scale);
        for(int p = 5; p <= 50  ; p += 5)
        {
            ArrayList<Long> seed_times = new ArrayList<>();
            ArrayList<Long> grow_times = new ArrayList<>();
            ArrayList<Long> enclaves_times = new ArrayList<>();
            ArrayList<Long> interregion_times = new ArrayList<>();
            ArrayList<Long> flow_times = new ArrayList<>();
            ArrayList<Long> local_times = new ArrayList<>();
            ArrayList<Long> before_heteros = new ArrayList<>();
            ArrayList<Long> after_heteros = new ArrayList<>();
            ArrayList<Double> rates = new ArrayList<>();
            int solved_cases = 0;

            for(int i = 0 ; i < GSLO_iter ; i++)
            {
                GlobalSearch sol = new GlobalSearch(Area.area_list_copy(all_areas) , p , all_areas.size(), threshold , false);
                seed_times.add(sol.getSeed_time());
                grow_times.add(sol.getRegion_growth_time());
                enclaves_times.add(sol.getEnclaves_assign_time());

                if(sol.isInterregion_flag())
                {
                    interregion_times.add(sol.getInterregion_update_time());
                }

                if(sol.isFlow_flag())
                {
                    flow_times.add(sol.getIndirect_flow_time());
                }

                if(sol.solved())
                {
                    solved_cases++;
                    long before_hetero = Region.get_all_region_hetero(sol.get_regions());
                    LocalOptimization up = new LocalOptimization(sol , all_areas.size() , 0.99 , sol.get_all_areas() , sol.get_regions() , threshold);
                    local_times.add(up.getTotal_time());
                    long after_hetero = up.getBest_hetero();
                    double improve = ((before_hetero - after_hetero) * 1.0) / before_hetero;
                    before_heteros.add(before_hetero);
                    after_heteros.add(after_hetero);
                    rates.add(improve);

                }



            }

            System.out.println("the current p is " + p + " the seed runtime is " + compute_long_ave(seed_times) + " region growth runtime " + compute_long_ave(grow_times) + " enclaves time " + compute_long_ave(enclaves_times));
            if(interregion_times.size() > 0)
            {
                System.out.println("num time entered interregion " + interregion_times.size() + "the interregional time is " + compute_long_ave(interregion_times));
            }

            if(flow_times.size() > 0)
            {
                System.out.println("num times entered flow " + flow_times.size() +   " the flow time is " + compute_long_ave(flow_times));
            }

            if(after_heteros.size() > 0)
            {
                System.out.println("avg before is " + compute_long_ave(before_heteros) + " ave after is " + compute_long_ave(after_heteros) + "rate is " + compute_double_ave(rates) + " local time is " + compute_long_ave(local_times));
            }

            System.out.println(" the solved cases is " + solved_cases);

        }
    }




    public static long compute_long_ave(ArrayList<Long> long_arr)
    {
        if(long_arr.size() == 0)
        {
            return -1;
        }
        long total = 0;
        for(long l : long_arr)
        {
            total += l;
        }

        return total / long_arr.size();
    }

    public static double compute_double_ave(ArrayList<Double> double_arr)
    {
        if(double_arr.size() == 0)
        {
            return  -1;
        }
        double total = 0;
        for(double d : double_arr)
        {
            total += d;
        }

        return total / double_arr.size();
    }


















}
