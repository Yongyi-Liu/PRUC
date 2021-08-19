package baseline.greedy;

import util.Area;

import util.Region;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;

/**
 * This class describes the hill climbing searching for the greedy baseline algorithm
 */
public class GreedySearch {
    private ArrayList<Area> all_areas;
    private Region[] regions;
    private long threshold;
    private long total_time;
    private long after_hetero;
    private SimpleGreedy sol;
    private long best_hetero;
    private int iter;
    private long start;

    /**
     *
     * @param sol The initial partition from the Simple Greedy
     * @param all_areas the input set of areas
     * @param regions the constructed regions
     * @param threshold the value of the user-defined constraint
     */
    public GreedySearch(SimpleGreedy sol , ArrayList<Area> all_areas , Region[] regions , long threshold)
    {
        long start = System.currentTimeMillis();
        this.start = start;
        this.sol = sol;
        this.all_areas = all_areas;
        this.regions = regions;
        this.threshold = threshold;
        this.after_hetero = Region.get_all_region_hetero(sol.getRegions());
        this.best_hetero = this.after_hetero;
        heuristics();
        this.total_time = System.currentTimeMillis() - start;
    }

    /**
     * The heuristics searches the move that brings the best heterogeneity improvement on the partition in each iteration.
     * When the number of iterations reaches the predefined value, it terminates.
     */
    public void heuristics()
    {
        if(!sol.is_solved())
        {
            this.total_time = 0;
            return;
        }
        while(true)
        {
            Move move = find_best_move();
            Region donor = move.getBelonging_region();
            Region receiver = move.getNew_region();
            Area area = move.getArea();
            long absolute_decre = donor.compute_hetero_decre(area) - receiver.compute_hetero_incre(area);
            donor.remove_area_in_region(area);
            receiver.add_area_to_region(area);
            this.after_hetero -= absolute_decre;
            iter++;
            if(this.after_hetero < this.best_hetero)
            {
                this.best_hetero = this.after_hetero;
            }


            if(iter == 200)
            {
                return;
            }


        }

    }

    /**
     * This method identifies the move that has the most reduction on the heterogeneity
     * @return the best move
     */
    public Move find_best_move(){
        ArrayList<Move> moves = new ArrayList<>();
        for(Area area : all_areas)
        {
            Region r = regions[area.get_associated_region_index()];
            ArrayList<Region> neigh_regions = find_neighbor_region(area);
            if(neigh_regions.size() > 0 && !r.area_disconect_region(area) && (r.get_region_extensive_attr() - area.get_extensive_attr() >= threshold))
            {
                for(Region neigh_r : neigh_regions)
                {
                    long hetero_decre = r.compute_hetero_decre(area);
                    long heter_incre = neigh_r.compute_hetero_incre(area);
                    moves.add(new Move(area , r , neigh_r , (hetero_decre - heter_incre)));
                }
            }
        }
        return Collections.max(moves, Comparator.comparingLong(Move::getHetero_decre));
    }

    public ArrayList<Region> find_neighbor_region(Area area)
    {
        ArrayList<Region> neighbor_regions = new ArrayList<>();
        Region r = regions[area.get_associated_region_index()];
        for(Area neigh_area : area.get_neigh_area(all_areas))
        {
            Region neigh_r = regions[neigh_area.get_associated_region_index()];
            if(r != neigh_r && !neighbor_regions.contains(neigh_r))
            {
                neighbor_regions.add(neigh_r);
            }
        }
        return neighbor_regions;
    }

    public long getBest_hetero()
    {
        return best_hetero;
    }


    public long getTotal_time()
    {
        return total_time;
    }
}
