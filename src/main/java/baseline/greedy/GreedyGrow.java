package baseline.greedy;

import GSLO.Seed;
import util.Area;
import util.Region;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;

/**
 * This class describes the greedy grow strategy for the greedy baseline for PRUC
 * The greedy grow approach grows the regions by adding the area that minimize the heterogeneity increase to the region
 */
public class GreedyGrow {
    private Seed seed;
    private Region[] regions;
    private long threshold;
    private ArrayList<Area> all_areas;
    private int p;

    public GreedyGrow(Seed seed , long threshold , ArrayList<Area> all_areas)
    {
        this.seed = seed;
        this.threshold = threshold;
        this.regions = new Region[seed.get_seed_size()];
        this.all_areas = all_areas;
        this.p = seed.get_seed_size();

    }

    /**
     * This method grows the regions sequentially
     * @return the set of grown regions
     */
    public Region[] grow_region_in_turn()
    {
        Comparator<Region> r_comparator = Comparator.comparingLong(r -> r.get_region_extensive_attr());

        for(int i = 0 ; i < regions.length ; i++)
        {
            Region r = new Region(i , seed.get_seeds().get(i), threshold , all_areas);
            regions[i] = r;
        }

        ArrayList<Region> growing_region = new ArrayList<>();
        Collections.addAll(growing_region, regions);

        growing_region.sort(r_comparator);

        while(growing_region.size() > 0)
        {
            Region region_to_grow = Collections.min(growing_region , r_comparator);
            grow(region_to_grow , growing_region);
        }

        return regions;
    }


    private void grow(Region r, ArrayList<Region> all_growing_regions)
    {
        if(r.get_region_extensive_attr() > threshold)
        {
            all_growing_regions.remove(r);
            return;
        }

        Area area_to_add = greedy_grow(r);

        if(area_to_add == null)
        {
            all_growing_regions.remove(r);
            return;
        }

        else
        {
            r.add_area_to_region(area_to_add);
        }

    }

    /**
     * This method grows a region by adding an area that minimizes the heterogeneity increase
     * @param r the growing region
     * @return the area that has the minimum heterogeneity increase
     */
    private Area greedy_grow(Region r)
    {
        ArrayList<Area> neighs = r.get_neigh_areas();
        long best_hetero_incre = Long.MAX_VALUE;
        Area best_area = null;
        for (Area current_area : neighs) {
            if (current_area.get_associated_region_index() != -1) {
                continue;
            }
            long hetero_incre = r.compute_hetero_incre(current_area);
            if (hetero_incre < best_hetero_incre) {
                best_hetero_incre = hetero_incre;
                best_area = current_area;
            }
        }
        return best_area;
    }



}
