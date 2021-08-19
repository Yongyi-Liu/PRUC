package util;

import org.geotools.data.DataStore;
import org.geotools.data.DataStoreFinder;
import org.geotools.data.FeatureSource;
import org.geotools.feature.FeatureCollection;
import org.geotools.feature.FeatureIterator;
import org.locationtech.jts.geom.Coordinate;
import org.locationtech.jts.geom.Geometry;
import org.opengis.feature.simple.SimpleFeature;
import org.opengis.feature.simple.SimpleFeatureType;
import org.opengis.filter.Filter;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Map;

/**
 * This class preprocess the data by reading the data from the shapefile, creating the area class and build the neighboring relations
 */
public class Preprocess {


    public static ArrayList<Area> GeoSetBuilder(String dataset) throws IOException {
        ArrayList<Area> areas = new ArrayList<>();
        FeatureCollection<SimpleFeatureType, SimpleFeature> collection = preprocess(dataset);
        ArrayList<Geometry> polygons = initial_construct(collection , areas, dataset);
        setNeighbors(polygons , areas);
        return areas;
    }

    private static FeatureCollection<SimpleFeatureType, SimpleFeature> preprocess(String dataset) throws IOException {

        File file = null;
        switch (dataset) {
            case "2k":
                file = new File("DataFile/2056dataset/merged.shp");
                break;
            case "diversity":
                file = new File("DataFile/diversity/2000data.shp");
                break;
            case "island":
                file = new File("DataFile/islanddata/WAandPENN.shp");
                break;
            case "5k":
                file = new File("DataFile/5K/5K.shp");
                break;
            case "10k":
                file = new File("DataFile/10K/10K.shp");
                break;
            case "20k":
                file = new File("DataFile/20K/20K.shp");
                break;
            case "30k":
                file = new File("DataFile/30K/30K.shp");
                break;
            case "40k":
                file = new File("DataFile/40K/40K.shp");
                break;
            case "50k":
                file = new File("DataFile/50K/50K.shp");
                break;
            case "60k":
                file = new File("DataFile/60K/60K.shp");
                break;
            case "70k":
                file = new File("DataFile/70K/70K.shp");
                break;
            case "80k":
                file = new File("DataFile/80K/80K.shp");
                break;
        }
        //System.out.println(file.getTotalSpace());
        Map<String, Object> map = new HashMap<>();
        map.put("url", file.toURI().toURL());
        DataStore dataStore = DataStoreFinder.getDataStore(map);
        String typeName = dataStore.getTypeNames()[0];
        FeatureSource<SimpleFeatureType, SimpleFeature> source =
                dataStore.getFeatureSource(typeName);
        Filter filter = Filter.INCLUDE;
        dataStore.dispose();
        return source.getFeatures(filter);
    }

    private static ArrayList<Geometry> initial_construct(FeatureCollection<SimpleFeatureType, SimpleFeature> collection , ArrayList<Area> areas, String dataset)
    {
        ArrayList<Geometry> polygons = new ArrayList<>();
        int geo_index = 0;
        try (FeatureIterator<SimpleFeature> features = collection.features()) {
            while (features.hasNext()) {
                SimpleFeature feature = features.next();
                long extensive_attr ;
                long internal_attr;

                if(dataset.equals("2k"))
                {
                    extensive_attr = Long.parseLong((feature.getAttribute("aland").toString()));
                    internal_attr  = Long.parseLong(feature.getAttribute("awater").toString());
                }

                else if(dataset.equals("diversity"))
                {
                    extensive_attr = (long)Double.parseDouble((feature.getAttribute("cty_pop200").toString()));
                    internal_attr = (long)(1000 * Double.parseDouble((feature.getAttribute("ratio").toString())));

                }

                else
                {
                    extensive_attr = Long.parseLong((feature.getAttribute("ALAND").toString()));
                    internal_attr  = Long.parseLong(feature.getAttribute("AWATER").toString());
                }

                Geometry polygon = (Geometry) feature.getDefaultGeometry();
                polygons.add(polygon);
                Coordinate[] coor = polygon.getCoordinates();
                Area newArea = new Area(geo_index , internal_attr , extensive_attr , coor);
                geo_index ++;
                areas.add(newArea);
            }
        }

        return polygons;

    }

    private static void setNeighbors(ArrayList<Geometry> polygons , ArrayList<Area> areas)
    {
        for (int i = 0; i < polygons.size(); i++) {

            for (int j = i + 1; j < polygons.size(); j++) {

                if (polygons.get(i).intersects(polygons.get(j))) {

                    Geometry intersection = polygons.get(i).intersection(polygons.get(j));

                    if (intersection.getGeometryType() != "Point") {

                        areas.get(i).add_neighbor(j);
                        areas.get(j).add_neighbor(i);


                    }
                }
            }
        }



    }
}
