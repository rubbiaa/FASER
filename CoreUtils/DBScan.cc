#include <iostream>
#include <vector>
#include <cmath>
#include <set>

#include "DBScan.hh"

// Density-Based Spatial Clustering of Applications with Noise (DBSCAN)

double DBScan::distance(const Point& p1, const Point& p2) {
    return std::sqrt(std::pow(p1.x - p2.x, 2) + std::pow(p1.y - p2.y, 2));
}

std::vector<int> DBScan::regionQuery(const std::vector<Point>& points, int index, double eps) {
    std::vector<int> neighbors;
    for (int i = 0; i < points.size(); ++i) {
        if (distance(points[index], points[i]) <= eps) {
            neighbors.push_back(i);
        }
    }
    return neighbors;
}

void DBScan::expandCluster(std::vector<Point>& points, int index, std::vector<int>& neighbors, int clusterID, double eps, int minPts) {
    points[index].clusterID = clusterID;
    int i = 0;
    while (i < neighbors.size()) {
        int neighborIndex = neighbors[i];
        if (!points[neighborIndex].visited) {
            points[neighborIndex].visited = true;
            std::vector<int> newNeighbors = regionQuery(points, neighborIndex, eps);
            if (newNeighbors.size() >= minPts) {
                neighbors.insert(neighbors.end(), newNeighbors.begin(), newNeighbors.end());
            }
        }
        if (points[neighborIndex].clusterID == -1) {
            points[neighborIndex].clusterID = clusterID;
        }
        ++i;
    }
}

void DBScan::scan(std::vector<Point>& points, double eps, int minPts) {
    int clusterID = 0;
    for (int i = 0; i < points.size(); ++i) {
        if (points[i].visited) continue;

        points[i].visited = true;
        std::vector<int> neighbors = regionQuery(points, i, eps);

        if (neighbors.size() < minPts) {
            points[i].clusterID = 0; // Mark as noise (could be any special ID, like 0)
        } else {
            ++clusterID;
            expandCluster(points, i, neighbors, clusterID, eps, minPts);
        }
    }
}


//////////////////////
// ---------- 3D DBSCAN ----------


double DBScan::distance3D(const Point3D& p1, const Point3D& p2) {
  return std::sqrt(std::pow(p1.x - p2.x, 2) + std::pow(p1.y - p2.y, 2) + std::pow(p1.z - p2.z, 2));
}


std::vector<int> DBScan::regionQuery3D(const std::vector<Point3D>& points, int index, double eps) {
  std::vector<int> neighbors;
  for (int i = 0; i < points.size(); ++i) {
    if (distance3D(points[index], points[i]) <= eps) {
      neighbors.push_back(i);
    }
  }
  return neighbors;
}

void DBScan::expandCluster3D(std::vector<Point3D>& points, int index, std::vector<int>& neighbors, int clusterID, double eps, int minPts) {
  points[index].clusterID = clusterID;
  int i = 0;
  while (i < neighbors.size()) {
    int neighborIndex = neighbors[i];
    if (!points[neighborIndex].visited) {
      points[neighborIndex].visited = true;
      std::vector<int> newNeighbors = regionQuery3D(points, neighborIndex, eps);
      if (newNeighbors.size() >= minPts) {
        neighbors.insert(neighbors.end(), newNeighbors.begin(), newNeighbors.end());
      }
    }
    if (points[neighborIndex].clusterID == -1) {
      points[neighborIndex].clusterID = clusterID;
    }
    ++i;
  }
}

void DBScan::scan3D(std::vector<Point3D>& points, double eps, int minPts) {
  int clusterID = 0;
  for (int i = 0; i < points.size(); ++i) {
    if (points[i].visited) continue;

    points[i].visited = true;
    std::vector<int> neighbors = regionQuery3D(points, i, eps);

    if (neighbors.size() < minPts) {
      points[i].clusterID = 0;
    } else {
      ++clusterID;
      expandCluster3D(points, i, neighbors, clusterID, eps, minPts);
    }
  }
}


#if 0
int main() {
    // Example points
    std::vector<Point> points = {
        {1.0, 2.0}, {2.0, 2.0}, {2.0, 3.0},
        {8.0, 7.0}, {8.0, 8.0}, {25.0, 80.0}
    };

    double eps = 2.0; // The maximum distance between two points to be considered as in the same neighborhood
    int minPts = 2;   // Minimum number of points to form a cluster

    dbscan(points, eps, minPts);

    // Output the cluster assignment
    for (const auto& point : points) {
        std::cout << "Point (" << point.x << ", " << point.y << ") -> Cluster " << point.clusterID << std::endl;
    }

    return 0;
}

#endif
