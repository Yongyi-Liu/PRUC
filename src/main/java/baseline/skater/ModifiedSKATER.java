package baseline.skater;

import util.Area;
import util.Region;

import java.util.ArrayList;
import java.util.Arrays;

/**
 * This class implements the Modified SKATER as one of our baseline algorithms
 * The original SKATER is presented in R. M. AssunÇão , M. C. Neves , G. Câmara & C. Da Costa Freitas (2006) Efficient regionalization techniques for socio‐economic geographical units using minimum spanning trees, International Journal of Geographical Information Science, 20:7, 797-811
 */
public class ModifiedSKATER {

    private int sc;
    private ArrayList<Area> all_areas;
    private int p;
    private long threshold;
    private long runtime;
    private ArrayList<Tree> trees;


    /**
     *
     * @param areas The input areas
     * @param sc the parameter that controls the stopping condition in the original SKATER paper
     * @param p the number of predefined regions
     * @param threshold the value on the user-defined constraint, when threshold is set to 0, the modified SKATER becomes SKATER
     * @param skatercon whether to run SKATER alone or run the SKATER for the input of SKATERCON
     */
    public ModifiedSKATER(ArrayList<Area> areas, int sc, int p, long threshold, boolean skatercon)
    {
        this.all_areas = areas;
        this.threshold = threshold;
        this.sc = sc;
        this.p = p;

        if(!skatercon)
        {
            long runtime_start = System.currentTimeMillis();
            Tree initial_tree = compute_MST_prim();
            trees = tree_partitioning(initial_tree);
            long runtime_end = System.currentTimeMillis();
            runtime = (runtime_end - runtime_start);
        }

        else
        {

            Tree initial_tree = new Tree(areas);
            trees = tree_partitioning(initial_tree);

        }

    }





    //we use prim's algorithm with priority queue implementation to compute the minimum spanning tree
    public Tree compute_MST_prim()
    {
        ArrayList<Edge> added_edges = new ArrayList<>();
        Area initial_area = all_areas.get(all_areas.size() - 1);
        int[] previous_index = new int[all_areas.size() - 1];
        Arrays.fill(previous_index , initial_area.get_geo_index());
        IndexMinPQ imp = new IndexMinPQ(all_areas.size() - 1);

        for(int i = 0 ; i < all_areas.size() - 1 ; i++)
        {
            if(initial_area.get_neigh_area_index().contains(i))
            {
                imp.insert(i+1 , initial_area.compute_hetero(all_areas.get(i)));
            }
            else
            {
                imp.insert(i+1 , Long.MAX_VALUE);
            }
        }

        for(int i = 0 ; i < all_areas.size() - 1 ; i++)
        {
            int index = imp.delMin() - 1;
            int previous = previous_index[index];
            added_edges.add(new Edge(index , previous));
            Area current_area = all_areas.get(index);
            for(int neighbor_index : current_area.get_neigh_area_index())
            {
                if(neighbor_index == initial_area.get_geo_index())
                {
                    continue;
                }

                long weight = current_area.compute_hetero(all_areas.get(neighbor_index));

                if(weight <  imp.keyOf(neighbor_index+1))
                {
                    imp.change(neighbor_index+1 , weight);
                    previous_index[neighbor_index] = index;
                }
            }
        }

        for(Area g : all_areas)
        {
            g.initialize_neighbor();
        }

        //reconstruct the neighborhood relationship using the neighborhood relationship stored in the edge
        for(Edge e : added_edges)
        {
            int from = e.get_from();
            int to = e.get_to();
            all_areas.get(from).add_neighbor(to);
            all_areas.get(to).add_neighbor(from);
        }

        return new Tree(all_areas);
    }


    /**
     * This method partitions the MST into p sub-trees, corresponding to the p regions
     * The f1 and f2 value are from the original SKATER paper
     * f1 is the heterogeneity of a tree and f2 is the heterogeneity different of two trees, which prevents generating two subtrees
     * that are very uneven from size
     * @param initial_tree
     * @return the p trees
     */
    private ArrayList<Tree> tree_partitioning(Tree initial_tree)
    {
        ArrayList<Tree> trees = new ArrayList<>();
        trees.add(initial_tree);

        int current_tree_num = 1;

        while(current_tree_num < p)
        {

            long best_f1 = 0;
            Edge best_edge = null;
            Tree t_to_partition = null;
            ArrayList<Area> best_t1 = null;
            ArrayList<Area> best_t2 = null;


            for(Tree t : trees)
            {
                if(t.getVertices().size() > 1)
                {
                    Object[] best_split_results = find_best_split(t);
                    long f1 = (long)best_split_results[0];
                    Edge edge = (Edge)best_split_results[1];

                    if(edge != null)
                    {
                        if(f1 > best_f1)
                        {
                            best_f1 = f1;
                            best_edge = edge;
                            best_t1 = (ArrayList<Area>)best_split_results[2];
                            best_t2 = (ArrayList<Area>)best_split_results[3];
                            t_to_partition = t;
                        }
                    }
                }
            }


            //if a feasible cut is not found among all the current existing trees, then the method returns null
            if(best_edge == null)
            {
                return null;
            }

            else
            {
                trees.remove(t_to_partition);
                Tree t1 = new Tree(best_t1);
                Tree t2 = new Tree(best_t2);
                trees.add(t1);
                trees.add(t2);
                int from = best_edge.get_from();
                int to = best_edge.get_to();
                all_areas.get(from).get_neigh_area_index().remove(Integer.valueOf(to));
                all_areas.get(to).get_neigh_area_index().remove(Integer.valueOf(from));
            }

            current_tree_num ++;

        }

        return trees;

    }

    /**
     * This method finds the best split, i.e. the split that results in a total of minimum heterogeniety
     * @param t the input tree
     * @return the best cut, and the split two sub-trees
     */
    private Object[] find_best_split (Tree t)
    {
        Area center_vertex = find_central_vertex(t);

        ArrayList<Area> visited = new ArrayList<>();

        Edge optimal_edge_split = null;
        long best_f1 = 0;
        int no_improve = 0;
        ArrayList<Area> best_t1 = null;
        ArrayList<Area> best_t2 = null;

        Area visiting_area = center_vertex;

        while(true)
        {

            visited.add(visiting_area);

            long best_neighbor_f2 = 0;
            Area next_to_expand = null;
            for(Area neigh : visiting_area.get_neigh_area(all_areas))
            {
                //make sure that we do not revisit a visited area
                if(visited.contains(neigh))
                {
                    continue;
                }
                Object[] split_results = split_by_edge(t , new Edge(visiting_area.get_geo_index() , neigh.get_geo_index()));
                boolean threshold_label = (boolean)split_results[0];
                long f1 = (long)split_results[1];
                long f2 = (long)split_results[2];
                ArrayList<Area> t1 = (ArrayList<Area>)split_results[3];
                ArrayList<Area> t2 = (ArrayList<Area>)split_results[4];

                if(f2 > best_neighbor_f2)
                {
                    best_neighbor_f2 = f2;
                    next_to_expand = neigh;
                }

                if(f1 > best_f1 && threshold_label)
                {
                    best_f1 = f1;
                    no_improve = 0;
                    optimal_edge_split = new Edge(visiting_area.get_geo_index() , neigh.get_geo_index());
                    best_t1 = t1;
                    best_t2 = t2;
                }
                else
                {
                    no_improve ++;
                }

                if(no_improve == sc)
                {
                    return new Object[]{best_f1 , optimal_edge_split , best_t1 , best_t2};
                }

            }
            if(next_to_expand == null)
            {
                return new Object[]{best_f1 , optimal_edge_split , best_t1 , best_t2};
            }

            visiting_area = next_to_expand;
        }
    }

    /**
     * this methonds evaluate the f1 and f2 values for a given tree when split by an edge
     * @param t The tree
     * @param e The edge to be removed from the tree
     * @return the two trees and the corresponding f1 f2 value
     */
    private Object[] split_by_edge(Tree t , Edge e)
    {
        ArrayList<Area> tree_vertices = t.getVertices();
        boolean[] visited = new boolean[tree_vertices.size()];
        DFS(tree_vertices.get(0) , tree_vertices , visited , e.get_from() , e.get_to());
        ArrayList<Area> t_1 = new ArrayList<>();
        ArrayList<Area> t_2 = new ArrayList<>();
        for(int i = 0 ; i < tree_vertices.size() ; i++)
        {
            if(visited[i])
            {
                t_1.add(tree_vertices.get(i));
            }
            else
            {
                t_2.add(tree_vertices.get(i));
            }
        }

        boolean threshold_label = true;

        if(compute_total_extensive(t_1) < threshold || compute_total_extensive(t_2) < threshold)
        {
            threshold_label = false;
        }

        long t1_hetero = compute_hetero(t_1);
        long t2_hetero = compute_hetero(t_2);

        long f1 = t.getTree_hetero()  - (t1_hetero + t2_hetero);
        long f2 = Math.max((t.getTree_hetero() - t1_hetero) , (t.getTree_hetero() - t2_hetero));

        return new Object[]{threshold_label , f1 , f2 , t_1 , t_2};
    }


    /**
     * This method finds the central vertex of a tree, i.e. the vertex that splits the tree to the closese two half when removed
     * @param t the tree to evaluate
     * @return the central area
     */
    private Area find_central_vertex(Tree t)
    {
        ArrayList<Area> vertices = t.getVertices();

        int half_size = vertices.size() / 2;

        ArrayList<Edge> visited_edges = new ArrayList<>();
        int best_diff = Integer.MAX_VALUE;
        Area best_split_vertex = null;

        for(Area v : vertices)
        {
            for(Area w : v.get_neigh_area(all_areas))
            {
                if(visited_edges.contains(new Edge(v.get_geo_index() ,w.get_geo_index())))
                {
                    continue;
                }
                visited_edges.add(new Edge(v.get_geo_index() , w.get_geo_index()));
                boolean[] visited = new boolean[vertices.size()];
                Area first_visit = vertices.get(0);
                DFS(first_visit , vertices , visited , v.get_geo_index(), w.get_geo_index());
                int count = 0;
                for(boolean visited_label : visited)
                {
                    if(visited_label)
                    {
                        count++;
                    }
                }

                if(Math.abs(half_size - count) < best_diff)
                {
                    best_diff = Math.abs(half_size - count);
                    best_split_vertex = v;

                }
            }

        }

        return best_split_vertex;


    }

    private void DFS(Area area_visit , ArrayList<Area> vertices , boolean[] visited , int from, int to)
    {
        visited[vertices.indexOf(area_visit)] = true;
        for(Area neigh : area_visit.get_neigh_area(all_areas))
        {
            if(!visited[vertices.indexOf(neigh)])
            {
                if((area_visit.get_geo_index() == from && neigh.get_geo_index() == to) || (area_visit.get_geo_index() == to && neigh.get_geo_index() == from))
                {
                    continue;
                }
                DFS(neigh , vertices , visited , from , to);
            }
        }
    }

    private long compute_hetero(ArrayList<Area> areas)
    {
        long hetero = 0;
        for(int i = 0; i < areas.size() ; i++)
        {
            for(int j = i + 1; j < areas.size() ; j++)
            {
                hetero += Math.abs(areas.get(i).get_internal_attr() - areas.get(j).get_internal_attr());
            }
        }
        return hetero;
    }


    private long compute_total_extensive(ArrayList<Area> areas)
    {
        long total_extensive = 0;
        for(Area area : areas)
        {
            total_extensive += area.get_extensive_attr();
        }
        return total_extensive;
    }

    public long getRuntime() {
        return runtime;
    }

    public ArrayList<Tree> getTrees()
    {
        return trees;
    }

    public Region[] get_regions()
    {
        if(trees == null)
        {
            return null;
        }
        Region[] regions = new Region[p];
        for(int i = 0 ; i < p ; i ++)
        {
            Tree t = trees.get(i);
            regions[i] = new Region(t.getVertices() , threshold , t.getTree_hetero() , t.getTree_total_extensive());
        }
        return regions;
    }





    class Tree
    {
        ArrayList<Area> vertices;
        long tree_hetero;
        long tree_total_extensive;

        public Tree(ArrayList<Area> vertices)
        {
            this.vertices = vertices;
            this.tree_hetero = compute_hetero(vertices);
            this.tree_total_extensive = compute_total_extensive(vertices);
        }

        public ArrayList<Area> getVertices()
        {
            return vertices;
        }

        public long getTree_hetero()
        {
            return tree_hetero;
        }

        public long getTree_total_extensive()
        {
            return tree_total_extensive;
        }

    }

    class Edge
    {
        int from_index;
        int to_index;

        public Edge(int from_index, int to_index)
        {
            this.from_index = from_index;
            this.to_index = to_index;
        }

        public int get_from()
        {
            return from_index;
        }

        public int get_to()
        {
            return to_index;
        }

        @Override
        public boolean equals(Object o) {
            if (this == o) return true;
            if (!(o instanceof Edge)) return false;
            Edge edge = (Edge) o;
            return (from_index == edge.from_index && to_index == edge.to_index) || (from_index == edge.to_index && to_index == edge.from_index);

        }


    }


    /**
     * The supporting datastructure to compute the MST
     */
    class IndexMinPQ {
        private int N;
        private int[] pq;
        private int[] qp;
        private long[] keys;

        public IndexMinPQ(int maxN) {
            keys =  new long[maxN + 1];
            pq = new int[maxN + 1];
            qp = new int[maxN + 1];
            Arrays.fill(qp, -1);
        }


        private boolean greater(int i, int j) {
            return keys[pq[i]] -keys[pq[j]] > 0;
        }

        public boolean isEmpty() {
            return N == 0;
        }

        public int size() {
            return N;
        }

        public boolean contains(int k) {
            return qp[k] != -1;
        }

        public void insert(int k, long key) {
            if (!contains(k)) {
                N++;
                pq[N] = k;
                qp[k] = N;
                keys[k] = key;
                swim(N);
            }
        }


        public void change(int k, long key) {
            keys[k] = key;
            swim(qp[k]);
            sink(qp[k]);
        }


        public long min() {
            return keys[pq[1]];
        }



        public int delMin() {
            int indexOfMin = pq[1];
            swap(1, N--);
            sink(1);
            keys[indexOfMin] = 0;

            qp[indexOfMin] = -1;

            return indexOfMin;
        }


        public long keyOf(int k) {
            if (contains(k)) {
                return keys[k];
            }
            return 0;
        }

        private void swap(int i, int j) {
            int temp = pq[i];
            pq[i] = pq[j];
            pq[j] = temp;
            qp[pq[i]] = i;
            qp[pq[j]] = j;
        }

        private void swim(int k) {
            while (k > 1 && greater(k / 2, k)) {
                swap(k / 2, k);
                k = k / 2;
            }
        }

        private void sink(int k) {
            while (2 * k <= N) {
                int j = 2 * k;
                if (j < N && greater(j, j + 1)) {
                    j++;
                }
                if (!greater(k, j)) {
                    break;
                }
                swap(k, j);
                k = j;
            }
        }
    }





}
