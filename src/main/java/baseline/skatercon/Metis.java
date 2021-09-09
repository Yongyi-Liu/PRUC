package baseline.skatercon;

import util.Area;
import util.Region;

import java.util.*;

/**
 * The graph partitioning algorithm to partition the similarity graph in SKATERCON.
 * The implementaion is based on Karypis, George, and Vipin Kumar. "A fast and high quality multilevel scheme for partitioning irregular graphs." SIAM Journal on scientific Computing 20, no. 1 (1998): 359-392.
 */
public class Metis {

    private ArrayList<Area> all_areas;
    private ArrayList<Region[]> regions_set;
    private int coarsen_threshold;
    private int p;
    private long threshold;

    public Metis(ArrayList<Area> all_areas_clean , ArrayList<Region[]> regions_set , int p , int coarsen_threshold , long threshold) throws CloneNotSupportedException {
        this.all_areas = Area.area_list_copy(all_areas_clean);
        this.regions_set = regions_set;
        this.coarsen_threshold = coarsen_threshold;
        this.p = p;
        this.threshold = threshold;
    }



    public Region[] Metis_start()
    {

        //if the threshold on the user-defined constraint is not 0, then the subgraph with the largest extensive attribute is given priority to be partitioned
        if(threshold > 0)
        {
            ArrayList<Graph> all_graphs = new ArrayList<>();
            Graph initial_graph = construct_similarity_graph(regions_set);
            ArrayList<Vertex> initial_vertices = initial_graph.getVertices();
            for(int i = 0 ; i < initial_vertices.size() ; i++)
            {
                Vertex v = initial_vertices.get(i);
                v.setGeo_index(i);
                v.setIndex(i);
            }


            all_graphs.add(initial_graph);

            while(all_graphs.size() < p)
            {
                Graph graph_to_split = select_largest_graph(all_graphs);
                all_graphs.remove(graph_to_split);
                ArrayList<Graph> graphs = coarsen(graph_to_split , coarsen_threshold);
                partition_min_graph(graphs.get(graphs.size() - 1));
                int flag = uncoarsen(graphs);
                if(flag == -1)
                {
                    return null;
                }

                Graph[] split_graph = split_graph(graph_to_split);
                if(split_graph == null)
                {
                    return null;
                }
                Graph left_graph = split_graph[0];
                Graph right_graph = split_graph[1];

                all_graphs.add(left_graph);
                all_graphs.add(right_graph);

            }


            Region[] regions = new Region[all_graphs.size()];

            int i = 0;
            while(i < p)
            {
                Graph g = all_graphs.get(i);
                ArrayList<Vertex> vertices = g.getVertices();
                ArrayList<Area> areas_in_region = new ArrayList<>();
                for(Vertex v : vertices)
                {
                    areas_in_region.add(all_areas.get(v.geo_index));
                }
                Region r = new Region(areas_in_region , threshold , 0 , 0);
                regions[i] = r;
                i++;
            }


            return regions;
        }


        else
        {
            Queue<Graph> all_graphs = new LinkedList<>();
            Graph initial_graph = construct_similarity_graph(regions_set);

            ArrayList<Vertex> initial_vertices = initial_graph.getVertices();
            for(int i = 0 ; i < initial_vertices.size() ; i++)
            {
                Vertex v = initial_vertices.get(i);
                v.setGeo_index(i);
                v.setIndex(i);
            }


            all_graphs.add(initial_graph);

            while(all_graphs.size() < p)
            {
                Graph graph_to_split = all_graphs.remove();
                ArrayList<Graph> graphs = coarsen(graph_to_split , coarsen_threshold);
                partition_min_graph(graphs.get(graphs.size() - 1));
                int flag = uncoarsen(graphs);
                if(flag == -1)
                {
                    return null;
                }

                Graph[] split_graph = split_graph(graph_to_split);
                if(split_graph == null)
                {
                    return null;
                }
                Graph left_graph = split_graph[0];
                Graph right_graph = split_graph[1];

                all_graphs.add(left_graph);
                all_graphs.add(right_graph);

            }


            Region[] regions = new Region[all_graphs.size()];

            int i = 0;
            while(i < p)
            {
                Graph g = all_graphs.remove();
                ArrayList<Vertex> vertices = g.getVertices();
                ArrayList<Area> areas_in_region = new ArrayList<>();
                for(Vertex v : vertices)
                {
                    areas_in_region.add(all_areas.get(v.geo_index));
                }
                Region r = new Region(areas_in_region , threshold, 0 , 0);
                regions[i] = r;
                i++;
            }


            return regions;
        }


    }







    private Graph construct_similarity_graph(ArrayList<Region[]> regions_set)
    {

        ArrayList<Vertex> vertices = new ArrayList<>();

        for(int i = 0 ; i < all_areas.size() ; i++)
        {
            Vertex v = new Vertex(1);
            v.setIndex(i);
            v.setGeo_index(i);
            vertices.add(v);
        }


        for(Region[] regions : regions_set)
        {
            for(Region r : regions)
            {
                ArrayList<Area> areas_in_r = r.get_areas_in_region();
                for(int i = 0 ; i < areas_in_r.size() ; i++)
                {
                    for(int j = i + 1 ; j < areas_in_r.size() ; j++)
                    {
                        Area g1 = all_areas.get(areas_in_r.get(i).get_geo_index());
                        Area g2 = all_areas.get(areas_in_r.get(j).get_geo_index());

                        if(g1.get_neigh_area_index().contains(g2.get_geo_index()))
                        {
                            int v1_index = g1.get_geo_index();
                            int v2_index = g2.get_geo_index();
                            Vertex v1 = vertices.get(v1_index);
                            Vertex v2 = vertices.get(v2_index);

                            if(v1.getNeighbor_index().contains(v2_index))
                            {
                                int v2_in_v1_pos = v1.getNeighbor_index().indexOf(v2_index);
                                v1.getEdge_weights().set(v2_in_v1_pos , v1.getEdge_weights().get(v2_in_v1_pos) + 1);
                                int v1_in_v2_pos = v2.getNeighbor_index().indexOf(v1_index);
                                v2.getEdge_weights().set(v1_in_v2_pos , v2.getEdge_weights().get(v1_in_v2_pos) + 1);
                            }

                            else
                            {
                                v1.getNeighbor_index().add(v2_index);
                                v2.getNeighbor_index().add(v1_index);
                                v1.getEdge_weights().add(1);
                                v2.getEdge_weights().add(1);
                            }

                        }
                    }
                }
            }


        }



        return new Graph(vertices);

    }


    private ArrayList<Graph> coarsen(Graph init_graph , int coarsen_threshold)
    {
        ArrayList<Graph> graphs = new ArrayList<>();
        Graph current_graph = init_graph;
        graphs.add(current_graph);

        for(int i = 0 ; i < current_graph.getVertices().size() ; i++)
        {
            current_graph.getVertices().get(i).setIndex(i);
        }




        while(current_graph.get_size() > coarsen_threshold)
        {
            Object[] matching_result = compute_maximal_matching(current_graph);
            int[] map = (int[])matching_result[0];
            int[] match = (int[])matching_result[1];
            current_graph.set_mapping_info(map , match);
            Graph contracted_graph = contract_graph(current_graph);
            graphs.add(contracted_graph);

            current_graph = contracted_graph;

        }

        return graphs;
    }

    private Graph contract_graph(Graph graph)
    {

        ArrayList<Vertex> previous_vertices = graph.getVertices();
        int[] map = graph.getMap();
        int[] match = graph.getMatch();

        int biggest_index = -1;
        for(int index : map)
        {
            if(index > biggest_index)
            {
                biggest_index = index;
            }
        }

        Vertex[] vertices_in_new_graph = new Vertex[biggest_index + 1];


        for(int v1 = 0 ; v1 < previous_vertices.size() ; v1++)
        {

            int u1 = map[v1];

            if(vertices_in_new_graph[u1] != null) //suggesting that u1 has already been constructed
            {
                continue;
            }

            ArrayList<Integer> u1_neighbors = new ArrayList<>();
            ArrayList<Integer> u1_neigh_weights = new ArrayList<>();

            int u1_weight = 0;
            ArrayList<Integer> v1_v2 = new ArrayList<>();
            v1_v2.add(v1);

            u1_weight += previous_vertices.get(v1).getVwget();


            int v2 = match[v1];

            if(v2 != -1)
            {
                u1_weight += previous_vertices.get(v2).getVwget();
                v1_v2.add(v2);
            }


            for(int v : v1_v2)
            {
                ArrayList<Integer> v_neighs = previous_vertices.get(v).getNeighbor_index();
                ArrayList<Integer> v_neigh_weights = previous_vertices.get(v).getEdge_weights();
                for(int i = 0 ; i < v_neighs.size() ; i++)
                {
                    int v_neigh = v_neighs.get(i);
                    if(map[v_neigh] == u1)
                    {
                        continue;
                    }

                    if(!u1_neighbors.contains(map[v_neigh]))
                    {
                        u1_neighbors.add(map[v_neigh]);
                        u1_neigh_weights.add(v_neigh_weights.get(i));
                    }

                    else
                    {

                        int index = (u1_neighbors.indexOf(map[v_neigh]));
                        {
                            u1_neigh_weights.set(index , u1_neigh_weights.get(index) + v_neigh_weights.get(i));
                        }
                    }
                }
            }

            Vertex new_v = new Vertex(u1 , u1_weight , u1_neighbors , u1_neigh_weights);
            vertices_in_new_graph[u1] = new_v;

        }

        ArrayList<Vertex> new_vertices = new ArrayList<>();
        Collections.addAll(new_vertices, vertices_in_new_graph);

        return new Graph(new_vertices);
    }


    private int[] choose_random_num(int n)
    {
        int len = n;
        int[] source = new int[len];
        for (int i = 0; i < len; i++)
        {
            source[i] = i;
        }
        int[] result = new int[n];
        Random rd = new Random(System.nanoTime());
        int index;
        for (int i = 0; i < result.length; i++)
        {
            index = Math.abs(rd.nextInt() % len--);
            result[i] = source[index];
            source[index] = source[len];
        }

        return result;
    }

    private Object[] compute_maximal_matching(Graph graph)
    {
        ArrayList<Vertex> vertices = graph.getVertices();
        int[] map = new int[vertices.size()];
        int[] match = new int[vertices.size()];
        int[] visiting_sequence = choose_random_num(vertices.size());

        Arrays.fill(map , -1);
        Arrays.fill(match , -1);

        int current_mapped_index = 0;

        for(int visiting_index : visiting_sequence)
        {
            if(map[visiting_index] != -1)
            {
                continue;
            }

            Vertex v = vertices.get(visiting_index);
            ArrayList<Integer> v_neighbor = v.getNeighbor_index();
            ArrayList<Integer> v_weights = v.getEdge_weights();

            int best_matching = -1;
            int max_weight = -1;

            for(int j = 0 ; j < v_neighbor.size() ; j++)
            {
                if(map[v_neighbor.get(j)] == -1)
                {
                    if(v_weights.get(j) > max_weight)
                    {
                        max_weight = v_weights.get(j);
                        best_matching = v_neighbor.get(j);
                    }
                }
            }





            map[visiting_index] = current_mapped_index;

            if(best_matching != -1)
            {
                map[best_matching] = current_mapped_index;
                match[visiting_index] = best_matching;
                match[best_matching] = visiting_index;
            }


            current_mapped_index ++;

        }

        return new Object[]{map , match};
    }



    public void partition_min_graph(Graph graph)
    {

        ArrayList<Vertex> vertices = graph.getVertices();


        boolean[] best_partition_label = null;
        int[] best_ID = new int[vertices.size()];
        int[] best_ED = new int[vertices.size()];

        int optimal_cross_weight = Integer.MAX_VALUE;

        for(int iter = 0 ; iter < 4 ; iter++)
        {
            Object[] partition_results = GGGP(vertices);
            int cross_edge_weight = (int)partition_results[0];
            if(cross_edge_weight < optimal_cross_weight)
            {
                optimal_cross_weight = cross_edge_weight;
                best_partition_label = (boolean[])partition_results[1];
                best_ED = (int[])partition_results[2];
                best_ID = (int[])partition_results[3];
            }
        }

        graph.setPartition_info(best_partition_label);
        graph.set_degree_info(best_ID ,best_ED);

    }




    private Object[] GGGP(ArrayList<Vertex> vertices)
    {
        boolean[] in_partition = new boolean[vertices.size()];
        ArrayList<Vertex> vertices_in_partition = new ArrayList<>();
        int current_total_weight = 0;
        int total_weight = 0;

        for(Vertex v : vertices)
        {
            total_weight += v.getVwget();
        }

        int half_total_weight = total_weight / 2;

        Vertex starting_v = vertices.get(new Random(System.nanoTime()).nextInt(vertices.size()));
        in_partition[starting_v.getIndex()] = true;
        vertices_in_partition.add(starting_v);

        current_total_weight += starting_v.getVwget();

        while(current_total_weight < half_total_weight)
        {
            ArrayList<Vertex> candidate_vertices = new ArrayList<>();
            ArrayList<Integer> potential_increase = new ArrayList<>();

            for(Vertex v : vertices_in_partition)
            {
                for(int v_candidate_index : v.getNeighbor_index())
                {
                    Vertex v_candidate = vertices.get(v_candidate_index);
                    if(!in_partition[v_candidate.getIndex()] && !candidate_vertices.contains(v_candidate))
                    {
                        candidate_vertices.add(v_candidate);

                        int weight_connecting_in = 0;
                        int weight_connecting_out = 0;

                        for(int i = 0 ; i < v_candidate.getNeighbor_index().size() ; i++)
                        {
                            if(in_partition[v_candidate.getNeighbor_index().get(i)])
                            {
                                weight_connecting_in += v_candidate.getEdge_weights().get(i);
                            }

                            else
                            {
                                weight_connecting_out += v_candidate.getEdge_weights().get(i);
                            }
                        }

                        int increase = weight_connecting_out - weight_connecting_in;

                        //in each iteration we add the vertex that mostly decreases the cross-partition weight
                        potential_increase.add(increase);

                    }
                }
            }

            Vertex best_vertex = null;

            if(candidate_vertices.size() == 0)
            {
                for(int i = 0 ; i < in_partition.length ; i++)
                {
                    if(!in_partition[i])
                    {
                        best_vertex = vertices.get(i);
                    }
                }
            }

            else
            {
                int min_index = potential_increase.indexOf(Collections.min(potential_increase));
                best_vertex = candidate_vertices.get(min_index);

            }

            in_partition[best_vertex.getIndex()] = true;
            current_total_weight += best_vertex.getVwget();
            vertices_in_partition.add(best_vertex);
        }

        int[] ED = new int[vertices.size()];
        int[] ID = new int[vertices.size()];

        for(int i = 0 ; i < vertices.size() ; i++)
        {
            int id = 0;
            int ed = 0;

            Vertex v = vertices.get(i);
            for(int j = 0 ; j < v.getNeighbor_index().size() ; j++)
            {
                int neigh = v.getNeighbor_index().get(j);
                if(in_partition[neigh] == in_partition[i])
                {
                    id += v.getEdge_weights().get(j);
                }

                else
                {
                    ed += v.getEdge_weights().get(j);
                }
            }

            ED[i] = ed;
            ID[i] = id;
        }

        int cross_edge_weight = 0;

        for(int ed : ED)
        {
            cross_edge_weight += ed;
        }

        cross_edge_weight = cross_edge_weight /2;

        return new Object[]{cross_edge_weight , in_partition, ED , ID};

    }


    private int BKL(Graph graph)
    {
        ArrayList<Vertex> vertices = graph.getVertices();
        int[] ED = graph.getED();
        int[] ID = graph.getID();
        boolean[] partition_label = graph.getPartition_label();

        int cross_edge_weight = 0;

        for(int ed : ED)
        {
            cross_edge_weight += ed;
        }

        cross_edge_weight = cross_edge_weight / 2;



        int best_cross_edge_weight = cross_edge_weight;
        int[] best_ED = ED.clone();
        int[] best_ID = ID.clone();
        boolean[] best_partition_label = partition_label.clone();

        int in_partition_weight = 0;
        int out_partition_weight = 0;
        Hashtable<Integer , Integer> in_partition_table = new Hashtable<>();
        Hashtable<Integer , Integer> out_partition_table = new Hashtable<>();

        for(int i = 0 ; i < vertices.size() ; i++)
        {
            Vertex v = vertices.get(i);
            if(partition_label[i])
            {
                in_partition_weight += v.getVwget();
                if(ED[i] > 0)
                {
                    in_partition_table.put(i , ED[i] - ID[i]);
                }

            }

            else
            {
                out_partition_weight += v.getVwget();
                if(ED[i] > 0)
                {
                    out_partition_table.put(i , ED[i] - ID[i]);
                }
            }
        }

        if(in_partition_table.isEmpty() || out_partition_table.isEmpty())
        {
            return -1;
        }


        boolean[] swapped = new boolean[vertices.size()];
        int no_improve = 0;
        while(no_improve < 50)
        {

            Vertex v;
            int index;
            if(in_partition_weight > out_partition_weight)
            {
                index = selecting_optimal_vertex(in_partition_table , swapped);
                if(index == -1)
                {
                    Arrays.fill(swapped , false);
                    continue;
                }
                v = vertices.get(index);
                in_partition_weight -= vertices.get(index).getVwget();
                out_partition_weight += vertices.get(index).getVwget();
                in_partition_table.remove(index);
            }

            else
            {
                index = selecting_optimal_vertex(out_partition_table , swapped);
                if(index == -1)
                {
                    Arrays.fill(swapped , false);
                    continue;
                }
                v = vertices.get(index);
                in_partition_weight += vertices.get(index).getVwget();
                out_partition_weight -= vertices.get(index).getVwget();
                out_partition_table.remove(index);
            }


            swapped[index] = true;
            partition_label[index] = (!partition_label[index]);


            int t = ED[index];
            ED[index] = ID[index];
            ID[index] = t;
            cross_edge_weight = cross_edge_weight + ED[index] - ID[index];

            if(partition_label[index])
            {
                in_partition_table.put(index , ED[index] - ID[index]);
            }
            else
            {
                out_partition_table.put(index, ED[index] - ID[index]);
            }



            for(int i = 0 ; i < v.getNeighbor_index().size() ; i++)
            {
                int neigh_v_index = v.getNeighbor_index().get(i);

                int edge_weight = v.getEdge_weights().get(i);
                if(partition_label[neigh_v_index] == partition_label[v.getIndex()])
                {
                    ED[neigh_v_index] -= edge_weight;
                    ID[neigh_v_index] += edge_weight;
                }
                else
                {
                    ED[neigh_v_index] += edge_weight;
                    ID[neigh_v_index] -= edge_weight;
                }

                if(partition_label[neigh_v_index])
                {
                    in_partition_table.remove(neigh_v_index);
                    if(ED[neigh_v_index] > 0)
                    {
                        in_partition_table.put(neigh_v_index , ED[neigh_v_index] - ID[neigh_v_index]);
                    }
                }

                else
                {
                    out_partition_table.remove(neigh_v_index);
                    if(ED[neigh_v_index] > 0)
                    {
                        out_partition_table.put(neigh_v_index , ED[neigh_v_index] - ID[neigh_v_index]);
                    }
                }

            }


            if(cross_edge_weight < best_cross_edge_weight)
            {
                no_improve = 0;
                best_cross_edge_weight = cross_edge_weight;
                best_ED = ED.clone();
                best_ID = ID.clone();
                best_partition_label = partition_label.clone();

            }

            else
            {
                no_improve ++;
            }

        }

        graph.set_degree_info(best_ID , best_ED);
        graph.setPartition_info(best_partition_label);



        return 0;


    }



    private int BKL1(Graph graph)
    {
        ArrayList<Vertex> vertices = graph.getVertices();
        int[] ED = graph.getED();
        int[] ID = graph.getID();
        boolean[] partition_label = graph.getPartition_label();

        int cross_edge_weight = 0;

        for(int ed : ED)
        {
            cross_edge_weight += ed;
        }

        cross_edge_weight = cross_edge_weight / 2;



        int best_cross_edge_weight = cross_edge_weight;
        int[] best_ED = ED.clone();
        int[] best_ID = ID.clone();
        boolean[] best_partition_label = partition_label.clone();

        int in_partition_weight = 0;
        int out_partition_weight = 0;
        Hashtable<Integer , Integer> in_partition_table = new Hashtable<>();
        Hashtable<Integer , Integer> out_partition_table = new Hashtable<>();

        for(int i = 0 ; i < vertices.size() ; i++)
        {
            Vertex v = vertices.get(i);
            if(partition_label[i])
            {
                in_partition_weight += v.getVwget();
                if(ED[i] > 0)
                {
                    in_partition_table.put(i , ED[i] - ID[i]);
                }

            }

            else
            {
                out_partition_weight += v.getVwget();
                if(ED[i] > 0)
                {
                    out_partition_table.put(i , ED[i] - ID[i]);
                }
            }
        }

        if(in_partition_table.isEmpty() || out_partition_table.isEmpty())
        {
            return -1;
        }

        boolean[] swapped = new boolean[vertices.size()];

        int no_improve = 0;

        while(no_improve < 50)
        {

            Vertex v;
            int index;
            if(in_partition_weight > out_partition_weight)
            {
                index = selecting_optimal_vertex(in_partition_table , swapped);
                if(index == -1)
                {
                    break;
                }
                v = vertices.get(index);
                in_partition_weight -= vertices.get(index).getVwget();
                out_partition_weight += vertices.get(index).getVwget();
                in_partition_table.remove(index);
            }

            else
            {
                index = selecting_optimal_vertex(out_partition_table , swapped);
                if(index == -1)
                {
                    break;
                }
                v = vertices.get(index);
                in_partition_weight += vertices.get(index).getVwget();
                out_partition_weight -= vertices.get(index).getVwget();
                out_partition_table.remove(index);
            }


            swapped[index] = true;
            partition_label[index] = (!partition_label[index]);


            int t = ED[index];
            ED[index] = ID[index];
            ID[index] = t;
            cross_edge_weight = cross_edge_weight + ED[index] - ID[index];

            if(partition_label[index])
            {
                in_partition_table.put(index , ED[index] - ID[index]);
            }
            else
            {
                out_partition_table.put(index, ED[index] - ID[index]);
            }

            for(int i = 0 ; i < v.getNeighbor_index().size() ; i++)
            {
                int neigh_v_index = v.getNeighbor_index().get(i);

                int edge_weight = v.getEdge_weights().get(i);
                if(partition_label[neigh_v_index] == partition_label[v.getIndex()])
                {
                    ED[neigh_v_index] -= edge_weight;
                    ID[neigh_v_index] += edge_weight;
                }
                else
                {
                    ED[neigh_v_index] += edge_weight;
                    ID[neigh_v_index] -= edge_weight;
                }

                if(partition_label[neigh_v_index])
                {
                    in_partition_table.remove(neigh_v_index);
                    if(ED[neigh_v_index] > 0)
                    {
                        in_partition_table.put(neigh_v_index , ED[neigh_v_index] - ID[neigh_v_index]);
                    }
                }

                else
                {
                    out_partition_table.remove(neigh_v_index);
                    if(ED[neigh_v_index] > 0)
                    {
                        out_partition_table.put(neigh_v_index , ED[neigh_v_index] - ID[neigh_v_index]);
                    }
                }

            }


            if(cross_edge_weight < best_cross_edge_weight)
            {
                no_improve = 0;
                best_cross_edge_weight = cross_edge_weight;
                best_ED = ED.clone();
                best_ID = ID.clone();
                best_partition_label = partition_label.clone();

            }

            else
            {
                no_improve ++;
            }

        }

        graph.set_degree_info(best_ID , best_ED);
        graph.setPartition_info(best_partition_label);



        return 0;
    }


    private int selecting_optimal_vertex(Hashtable<Integer , Integer> ht , boolean[] processed)
    {
        int optimal_index = -1;
        int optimal_gain = Integer.MIN_VALUE;
        for(int index :ht.keySet())
        {
            int gain = ht.get(index);
            if(gain > optimal_gain && !processed[index])
            {
                optimal_gain = gain;
                optimal_index = index;
            }
        }
        return optimal_index;
    }






    private int uncoarsen(ArrayList<Graph> graphs)
    {

        int graph_index = graphs.size() - 1;
        while(graph_index != 0)
        {
            Graph current_graph = graphs.get(graph_index);
            ArrayList<Vertex> current_vertices = current_graph.getVertices();
            boolean[] current_partition = current_graph.getPartition_label();


            int[] current_ED = current_graph.getED();
            int[] current_ID = current_graph.getID();

            Graph last_graph = graphs.get(graph_index - 1);
            int[] last_map = last_graph.getMap();
            int[] last_match = last_graph.getMatch();
            ArrayList<Vertex> last_vertices = last_graph.getVertices();
            int[] last_ED = new int[last_match.length];
            int[] last_ID = new int[last_match.length];
            boolean[] last_partition = new boolean[last_match.length];


            //this deals with the partition lable in the previous partition label;
            for(int i = 0 ; i < last_partition.length ; i++)
            {
                last_partition[i] = current_partition[last_map[i]];
            }



            for(int i = 0 ; i < last_ED.length ; i++)
            {
                int mapped = last_map[i];
                Vertex v = current_vertices.get(mapped);
                Vertex v1 = last_vertices.get(i);
                int matched_id = last_match[i];

                if(current_ED[mapped] == 0)
                {
                    last_ED[i] = 0;
                    int id = 0;
                    for(int weight : v1.getEdge_weights())
                    {
                        id += weight;
                    }
                    last_ID[i] = id;


                    if(matched_id != -1)
                    {
                        Vertex v2 = last_vertices.get(matched_id);
                        last_ED[matched_id] = 0;

                        int v2_id = 0;
                        for(int weight : v2.getEdge_weights())
                        {
                            v2_id += weight;
                        }
                        last_ID[matched_id] = v2_id;
                    }





                }

                else if(current_ID[mapped] == 0)
                {
                     if(matched_id == -1) //suggesting that this node is a single node;
                     {
                         last_ID[i] = 0;
                         int ed = 0;
                         for(int weight : v1.getEdge_weights())
                         {
                             ed += weight;
                         }

                         last_ED[i] = ed;

                     }

                     else //suggesting that this node is a multinode
                     {
                         Vertex v2 = last_vertices.get(matched_id);
                         int index = v1.getNeighbor_index().indexOf(matched_id);
                         int v1_v2_weight = v1.getEdge_weights().get(index);
                         last_ID[i] = v1_v2_weight;
                         last_ID[matched_id] = v1_v2_weight;
                         int v1_ed = 0;
                         for(int weight : v1.getEdge_weights())
                         {
                             v1_ed += weight;
                         }
                         v1_ed -= v1_v2_weight;
                         last_ED[i] = v1_ed;

                         int v2_ed = 0;
                         for(int weight : v2.getEdge_weights())
                         {
                             v2_ed += weight;
                         }
                         v2_ed -= v1_v2_weight;
                         last_ED[matched_id] = v2_ed;

                     }



                }


                else
                {
                    int v1_ed = 0;
                    int v1_id = 0;
                    for(int j = 0 ; j < v1.getNeighbor_index().size() ; j++)
                    {
                        int v1_neigh = v1.getNeighbor_index().get(j);
                        if(last_partition[i] != last_partition[v1_neigh])
                        {
                            v1_ed += v1.getEdge_weights().get(j);
                        }
                        else
                        {
                            v1_id += v1.getEdge_weights().get(j);
                        }
                    }

                    last_ED[i] = v1_ed;
                    last_ID[i] = v1_id;

                    if(matched_id != -1)
                    {
                        int index = v1.getNeighbor_index().indexOf(matched_id);
                        int v1_v2_weight = v1.getEdge_weights().get(index);

                        last_ED[matched_id] = current_ED[v.getIndex()] - last_ED[i];
                        last_ID[matched_id] = current_ID[v.getIndex()] - last_ID[i] + 2 * v1_v2_weight;

                    }


                }
            }


            last_graph.set_degree_info(last_ID , last_ED);
            last_graph.setPartition_info(last_partition);
            int flag;
            if(last_graph.get_size() < (int)(0.02 * all_areas.size()))
            {
                flag = BKL1(last_graph);
            }

            else
            {
                flag = BKL(last_graph);
            }

            if(flag == -1)
            {
                return -1;
            }

            graph_index --;


        }
        return 0;


    }


    //this method splits a graph into two
    private Graph[] split_graph(Graph graph)
    {
        boolean[] partition = graph.getPartition_label();
        ArrayList<Vertex> vertices = graph.getVertices();


        ArrayList<Vertex> in_partition = new ArrayList<>();
        ArrayList<Vertex> out_partition = new ArrayList<>();

        for(int i = 0 ; i < vertices.size() ; i++)
        {
            if(partition[i])
            {
                in_partition.add(vertices.get(i));
            }

            else
            {
                out_partition.add(vertices.get(i));
            }
        }



        for (int i = 0 ; i < vertices.size() ; i++) {
            Vertex v = vertices.get(i);

            if(partition[i])
            {
                int index_in_partition = in_partition.indexOf(v);
                v.setIndex(index_in_partition);
            }

            else
            {
                int index_out_partition = out_partition.indexOf(v);
                v.setIndex(index_out_partition);
            }

            ArrayList<Integer> updated_index = new ArrayList<>();
            ArrayList<Integer> updated_weight = new ArrayList<>();

            ArrayList<Integer> v_neighbor = v.getNeighbor_index();
            ArrayList<Integer> v_weights = v.getEdge_weights();

            for (int j = 0; j < v_neighbor.size(); j++) {
                int neigh = v_neighbor.get(j);
                if (partition[i] == partition[neigh]) {
                    if(partition[i])
                    {
                        int index = in_partition.indexOf(vertices.get(neigh));
                        updated_index.add(index);
                        updated_weight.add(v_weights.get(j));
                    }

                    if(!partition[i])
                    {
                        int index = out_partition.indexOf(vertices.get(neigh));
                        updated_index.add(index);
                        updated_weight.add(v_weights.get(j));
                    }
                }
            }

            v.set_updated_info(updated_index , updated_weight);
        }

        long left_total_extensive = 0L;

        for(Vertex v : in_partition)
        {
            left_total_extensive += all_areas.get(v.getGeo_index()).get_extensive_attr();
        }

        if(left_total_extensive < threshold)
        {
            return null;
        }

        long right_total_extensive = 0L;
        for(Vertex v : out_partition)
        {
            right_total_extensive += all_areas.get(v.getGeo_index()).get_extensive_attr();
        }
        if(right_total_extensive < threshold)
        {
            return null;
        }
        //System.out.println("the left is " + left_total_extensive + " the right is " + right_total_extensive + " the threshold is " + threshold);

        Graph split_left = new Graph(in_partition);
        Graph split_right = new Graph(out_partition);

        return new Graph[]{split_left , split_right};


    }


    private Graph select_largest_graph(ArrayList<Graph> graphs)
    {

        Graph maximum_graph = null;
        long max_extensive_attr = 0L;

        for(Graph g : graphs)
        {
            long total_extensive_attr = 0L;
            for(Vertex vertex : g.getVertices())
            {
                total_extensive_attr += all_areas.get(vertex.getGeo_index()).get_extensive_attr();
            }
            if(total_extensive_attr > max_extensive_attr)
            {
                maximum_graph = g;
                max_extensive_attr = total_extensive_attr;
            }
        }



        return maximum_graph;
    }








    class Vertex implements Cloneable
    {
        int geo_index;
        int vwget;
        int index;

        ArrayList<Integer> neighbor_index;
        ArrayList<Integer> edge_weights;

        public Vertex(int vwget)
        {
            this.vwget = vwget;
            neighbor_index = new ArrayList<>();
            edge_weights = new ArrayList<>();
        }

        public Vertex(int index , int vwget , ArrayList<Integer> neighbor_index , ArrayList<Integer> edge_weights)
        {
            this.index = index;
            this.vwget = vwget;
            this.neighbor_index =  neighbor_index;
            this.edge_weights = edge_weights;
        }

        public void set_updated_info(ArrayList<Integer> new_neighbors , ArrayList<Integer> new_weights)
        {
            this.neighbor_index = new_neighbors;
            this.edge_weights = new_weights;
        }

        public void setIndex(int index)
        {
            this.index = index;
        }

        public void setGeo_index(int geo_index)
        {
            this.geo_index = geo_index;
        }

        public int getVwget()
        {
            return vwget;
        }

        public int getIndex() { return index; }

        public int getGeo_index() { return geo_index; }

        public ArrayList<Integer> getNeighbor_index()
        {
            return neighbor_index;
        }

        public ArrayList<Integer> getEdge_weights()
        {
            return edge_weights;
        }

    }

    class Graph
    {
        ArrayList<Vertex> vertices;
        int[] map;
        int[] match;
        boolean[] partition_label;
        int[] ID;
        int[] ED;

        public Graph(ArrayList<Vertex> vertices)
        {
            this.vertices = vertices;

        }

        public void set_mapping_info(int[] map , int[] match)
        {
            this.map = map;
            this.match = match;
        }

        public void setPartition_info(boolean[] partition_label)
        {
            this.partition_label = partition_label;
        }

        public void set_degree_info(int[] ID , int[] ED)
        {
            this.ID = ID;
            this.ED = ED;
        }




        public ArrayList<Vertex> getVertices()
        {
            return vertices;
        }

        public int[] getMap()
        {
            return map;
        }

        public int[] getMatch()
        {
            return match;
        }

        public int get_size()
        {
            return vertices.size();
        }


        public boolean[] getPartition_label()
        {
            return partition_label;
        }

        public int[] getID()
        {
            return ID;
        }

        public int[] getED()
        {
            return ED;
        }

    }









}
