package GSLO;


import util.Area;

import java.util.ArrayList;
import java.util.Comparator;
import java.util.HashSet;
import java.util.Random;


/**
 * This class corresponds to Section 5.1.1 Seed Identification
 */

public class SeedIdentification {
    private ArrayList<Area> all_geoareas;
    private Seed best_seed;


    /**
     *
     * @param all_areas the input areas
     * @param seed_num the number of seeded area, which equals to the number of predefined regions, p
     * @param max_iter the maximum number of iteration in Seed Identification, when set to 0 means random seeding
     * @param kmeanspp whether or not to use the k-mean ++ seeding
     * @param detect_island whether or not to consider island areas in the dataset
     */
    public SeedIdentification(ArrayList<Area> all_areas , int seed_num , int max_iter , boolean kmeanspp , boolean detect_island) {
        this.all_geoareas = all_areas;
        if(!detect_island)
        {
            if(!kmeanspp)
            {
                this.best_seed = naive_seed_selection(all_areas , seed_num , max_iter);
            }

            else
            {
                this.best_seed = kmeanspp(all_areas , seed_num);
            }
        }

        else
        {
            this.best_seed = island_seeding(all_areas , seed_num , max_iter);
        }
    }


    /**
     *
     * @param all_areas the input areas
     * @param s_num the number of areas in the seed
     * @param maxiter the maximum number of iterations
     * @return The selected seed
     */
    public Seed island_seeding(ArrayList<Area> all_areas , int s_num , int maxiter)
    {
        ArrayList<ConnectedComponent> ccs = new ArrayList<>();
        Comparator<ConnectedComponent> cc_comparator = Comparator.comparingLong(ConnectedComponent::getTotal_ext);

        ArrayList<Area> unvisited = (ArrayList<Area>) (all_areas.clone());
        while(unvisited.size() > 0)
        {
            Area a = unvisited.get(new Random().nextInt(unvisited.size()));
            ArrayList<Area> visited = find_connected_component(a);
            long cc_ext = 0;
            for(Area cc_a : visited)
            {
                cc_ext += cc_a.get_extensive_attr();
            }
            ccs.add(new ConnectedComponent(visited, cc_ext));
            unvisited.removeAll(visited);
        }

        ccs.sort(cc_comparator);

        long total_ext = 0;
        for(ConnectedComponent cc : ccs)
        {
            total_ext += cc.getTotal_ext();
        }


        int remaining_seed = s_num;
        ArrayList<Seed> seeds_in_cc = new ArrayList<>();

        for(int i = 0 ; i < ccs.size() - 1 ; i++)
        {
            ConnectedComponent c = ccs.get(i);
            int seed_num = (int) ((c.getTotal_ext() * 1.0 / total_ext) * s_num);
            if(seed_num < 1)
            {
                seed_num = 1;
            }
            int max_iter = (int) ((c.getTotal_ext() * 1.0 / total_ext) * maxiter);
            Seed s = naive_seed_selection(c.getAreas_in_cc() , seed_num , max_iter);
            seeds_in_cc.add(s);
            remaining_seed -= seed_num;
        }

        ConnectedComponent c = ccs.get(ccs.size() - 1);
        int max_iter = (int) ((c.getTotal_ext() * 1.0 / total_ext) * maxiter);
        Seed s = naive_seed_selection(c.getAreas_in_cc() , remaining_seed ,max_iter);
        seeds_in_cc.add(s);

        ArrayList<Area> all_seed_areas = new ArrayList<>();
        for(Seed s_in_cc : seeds_in_cc)
        {
            all_seed_areas.addAll(s_in_cc.get_seeds());
        }
        return new Seed(all_seed_areas);
    }

    public ArrayList<Area> find_connected_component(Area area)
    {
        HashSet<Area> visited = new HashSet<>();
        DFS(area , visited);
        return new ArrayList<>(visited);
    }

    public void DFS(Area area , HashSet<Area> visited)
    {
        visited.add(area);
        for(Area neigh_area : area.get_neigh_area(all_geoareas))
        {
            if(!visited.contains(neigh_area))
            {
                DFS(neigh_area , visited);
            }
        }
    }


    public Seed naive_seed_selection(ArrayList<Area> all_areas , int s_num , int maxiter)
    {
        Seed seed = new Seed(all_areas , s_num);
        int iter_time = 0;
        while(iter_time < maxiter)
        {
            seed.random_replacement();
            iter_time ++;
        }
        return seed;
    }


    public Seed kmeanspp(ArrayList<Area> all_areas , int s_num)
    {
        Area init_area = all_areas.get(new Random().nextInt(all_areas.size()));
        ArrayList<Area> seeded_areas = new ArrayList<>();
        seeded_areas.add(init_area);
        while(seeded_areas.size() < s_num)
        {
            ArrayList<Area> unseeded_areas = new ArrayList<>();
            ArrayList<Double> dists = new ArrayList<>();
            double total_dist = 0;
            for(Area a : all_areas)
            {
                if(!seeded_areas.contains(a))
                {
                    double min_dist_to_seed = Long.MAX_VALUE;
                    for(Area seed_a :seeded_areas)
                    {
                        double dist = a.compute_dist(seed_a);
                        if(dist < min_dist_to_seed)
                        {
                            min_dist_to_seed = dist;
                        }
                    }
                    unseeded_areas.add(a);
                    dists.add(min_dist_to_seed);
                    total_dist += (min_dist_to_seed * min_dist_to_seed);
                }
            }

            double random_num = new Random().nextDouble();
            double accu = 0;
            Area selected_area = null;
            for(int i = 0 ; i < unseeded_areas.size() ; i++)
            {
                double d = dists.get(i);
                double prob = d * d / total_dist;
                if(random_num >= accu && random_num <= accu + prob)
                {
                    selected_area = unseeded_areas.get(i);
                    break;
                }
                accu += prob;
            }



            if(selected_area == null)
            {
                selected_area = seeded_areas.get(seeded_areas.size() - 1);
            }
            seeded_areas.add(selected_area);

        }

        return new Seed(seeded_areas);
    }

    public Seed getBest_seed() {
        return best_seed;
    }


    static class ConnectedComponent{
        ArrayList<Area> areas_in_cc;
        long total_ext;
        public ConnectedComponent(ArrayList<Area> areas , long total_e)
        {
            areas_in_cc = areas;
            total_ext = total_e;
        }

        public ArrayList<Area> getAreas_in_cc()
        {
            return areas_in_cc;
        }

        public long getTotal_ext()
        {
            return total_ext;
        }

    }
}

