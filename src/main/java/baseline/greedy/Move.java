package baseline.greedy;

import util.Area;

import util.Region;

public class Move {
    Area area;
    Region belonging_region;
    Region new_region;
    long hetero_decre;


    public Move(Area area , Region belonging_region , Region new_region , long hetero_decre)
    {
        this.area = area;
        this.belonging_region = belonging_region;
        this.new_region = new_region;
        this.hetero_decre = hetero_decre;
    }

    public Area getArea()
    {
        return area;
    }

    public Region getBelonging_region()
    {
        return belonging_region;
    }

    public Region getNew_region()
    {
        return new_region;
    }

    public long getHetero_decre() {return hetero_decre;}
}
