package baseline.skatercon;
import util.Area;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Random;

/**
 * This class constructs a random spanning tree, which is the prerequisite of skatercon
 * From the paper David Bruce Wilson. 1996. Generating random spanning trees more quickly than the cover time. In Proceedings ofthe twenty-eighth annual ACM symposium on Theory ofcomputing. ACM, 296–303.
 */
public class RandomSpanningTree {

    private ArrayList<Area> all_areas;

    public RandomSpanningTree(ArrayList<Area> all_areas)
    {
        this.all_areas = (ArrayList<Area>) all_areas.clone();
    }

    public ArrayList<Area> run_RST()
    {
        boolean[] in_tree = new boolean[all_areas.size()];
        int[] next = new int[all_areas.size()];
        Area root = all_areas.get(new Random().nextInt(all_areas.size()));
        int root_index = all_areas.indexOf(root);
        in_tree[root_index] = true;
        Arrays.fill(next , -1);

        for(int i = 0 ; i < all_areas.size() ; i++)
        {
            int u = i;
            while(!in_tree[u])
            {
                next[u] = random_successor(u);
                u = next[u];
            }
            u = i;
            while(!in_tree[u])
            {
                in_tree[u] = true;
                u = next[u];
            }
        }

        for(Area g : all_areas)
        {
            g.initialize_neighbor();
        }

        for(int i = 0 ; i < all_areas.size() ; i++)
        {
            int j = next[i];
            if(j == -1)
            {
                continue;
            }
            all_areas.get(i).add_neighbor(j);
            all_areas.get(j).add_neighbor(i);

        }

        int total_neigh = 0;
        for(Area area : all_areas)
        {
            total_neigh +=  area.get_neigh_area_index().size();
        }

        return all_areas;


    }




    private int random_successor(int u)
    {
        ArrayList<Integer> neigh_index =  all_areas.get(u).get_neigh_area_index();
        int rand_num = new Random(System.nanoTime()).nextInt(neigh_index.size());
        return neigh_index.get(rand_num);

        /*
        GeoArea g = all_areas.get(u);
        ArrayList<Integer> all_neighbor_index = g.get_neigh_area_index();
        ArrayList<Double> all_neigh_inverse_weight = new ArrayList<>();
        double total_inverse_weight = 0;
        for(int neigh_index : all_neighbor_index)
        {
            long diff = Math.abs(g.get_internal_attr() - all_areas.get(neigh_index).get_internal_attr());
            double diff_inverse;
            if(diff == 0)
            {
                diff_inverse = 1;
            }
            else
            {
                diff_inverse = (1.0 / diff);
            }
            total_inverse_weight += diff_inverse;
            all_neigh_inverse_weight.add(diff_inverse);
        }

        ArrayList<Double> prob_interval = new ArrayList<>();
        double probability_accu = 0;
        prob_interval.add(0.0);
        for(int i = 0 ; i < all_neigh_inverse_weight.size() - 1; i++)
        {
            probability_accu += (all_neigh_inverse_weight.get(i) / total_inverse_weight);
            prob_interval.add(probability_accu);
        }
        prob_interval.add(1.0);

        double random_num = new Random(System.nanoTime()).nextDouble();




        for(int i = 0 ; i < prob_interval.size() ; i++)
        {
            if(random_num >= prob_interval.get(i) && random_num <= prob_interval.get(i+1))
            {
                return all_neighbor_index.get(i);
            }
        }

        return  all_neighbor_index.get(all_neighbor_index.size() - 1);


         */

    }




}
