package baseline.greedy;

import GSLO.EnclavesAssignment;
import GSLO.Seed;
import GSLO.SeedIdentification;
import util.Area;
import util.Region;

import java.util.ArrayList;
import java.util.HashSet;
import java.util.Random;

/**
 * This class corresponds to the simple greedy baseline to solve PRUC in Section 7
 */
public class SimpleGreedy {
    private Region[] regions;
    private ArrayList<Area> all_areas;
    private Seed seed;
    private long total_running_time;
    private boolean solved;
    private boolean island_no_seed;

    /**
     *
     * @param all_areas The input areas
     * @param p The predefined number of regions
     * @param threshold The value on the user-defined constraint
     * @param detect_island Whether or not to consider island, just for testing purposes here
     */
    public SimpleGreedy(ArrayList<Area> all_areas, int p, long threshold , boolean detect_island) {
        long start = System.currentTimeMillis();

        this.all_areas = all_areas;
        for (Area all_area : all_areas) { all_area.set_centroid(); }

        //when max_iter is set to 0, this corresponds to random seeding
        this.seed = new SeedIdentification(all_areas , p , 0 , false , false).getBest_seed();

        if(!detect_island)
        {
            this.regions = new GreedyGrow(seed , threshold , all_areas).grow_region_in_turn();
            new EnclavesAssignment(all_areas , this.regions);
            long end = System.currentTimeMillis();
            total_running_time = end - start;
            this.solved = is_solved();
        }

        else //if island presents, we need to exam whether each connected component has at least one seeded area
        {
            island_no_seed = check_island_seed(this.seed , all_areas);
            if(island_no_seed)
            {
                long end = System.currentTimeMillis();
                total_running_time = end - start;
                this.solved = false;
            }

            else
            {
                this.regions = new GreedyGrow(seed , threshold , all_areas).grow_region_in_turn();
                new EnclavesAssignment(all_areas , this.regions);
                long end = System.currentTimeMillis();
                total_running_time = end - start;
                this.solved = is_solved();
            }


        }



    }

    /**
     * This method checks whether each connected component has at least one seeded area
     * @param s the seed obtained by random sampling
     * @param all_areas the input areas
     * @return whether or not each connected component has at least one seeded area
     */
    public boolean check_island_seed(Seed s , ArrayList<Area> all_areas)
    {
        ArrayList<Area> unvisited = (ArrayList<Area>) (all_areas.clone());

        while(unvisited.size() > 0)
        {
            Area a = unvisited.get(new Random().nextInt(unvisited.size()));
            ArrayList<Area> visited = find_connected_component(a);
            ArrayList<Area> seed_copy = (ArrayList<Area>) s.get_seeds().clone();
            seed_copy.removeAll(visited);
            if(seed_copy.size() == s.get_seeds().size())
            {
                return true;
            }


            unvisited.removeAll(visited);
        }

        return false;


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
        for(Area neigh_area : area.get_neigh_area(all_areas))
        {
            if(!visited.contains(neigh_area))
            {
                DFS(neigh_area , visited);
            }
        }
    }


    public boolean is_solved()
    {
        if(regions == null)
        {
            return false;
        }
        for(Region r : regions)
        {
            if(!r.is_region_complete())
            {
                return false;
            }
        }
        return true;
    }

    public boolean get_solved()
    {
        return solved;
    }

    public ArrayList<Area> getAll_areas()
    {
        return all_areas;
    }

    public Region[] getRegions(){
        return regions;
    }

    public long getTotal_running_time()
    {
        return total_running_time;
    }
}
