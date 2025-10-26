#ifndef _DBSCAN_
#define _DBSCAN_ 1

#include <TObject.h>


/// @brief Density-Based Spatial Clustering of Applications with Noise (DBSCAN)
class DBScan : public TObject {
public:

struct Point {
    long ID;
    double ehit;
    double x, y;
    bool visited = false;
    int clusterID = -1; // -1 means not yet assigned to any cluster
};

 // Added by Umut
struct Point3D {
    long ID;
    double ehit;
    double x, y, z;
    int clusterID = -1;
    bool visited = false;
  };

    DBScan() {};
    virtual ~DBScan() {};

    void scan(std::vector<Point>& points, double eps, int minPts);
  void scan3D(std::vector<Point3D>& points, double eps, int minPts);

private:

    double distance(const Point& p1, const Point& p2);
    std::vector<int> regionQuery(const std::vector<Point>& points, int index, double eps);
    void expandCluster(std::vector<Point>& points, int index, std::vector<int>& neighbors, int clusterID, double eps, int minPts);
 double distance3D(const Point3D& p1, const Point3D& p2);
  std::vector<int> regionQuery3D(const std::vector<Point3D>& points, int index, double eps);
  void expandCluster3D(std::vector<Point3D>& points, int index, std::vector<int>& neighbors, int clusterID, double eps, int minPts);


};

#endif
