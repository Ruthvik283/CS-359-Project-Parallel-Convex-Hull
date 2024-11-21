#include <iostream>
#include <vector>
#include <algorithm>
#include <utility>
#include <cmath>
#include <chrono>
#include <omp.h>

using namespace std;

using ll = long long;
using Point = pair<ll, ll>;  // Pair representing a point (x, y)
using VectorPoints = vector<Point>;

ll square(ll a) { return a * a; } 
ll norm(const Point &p) { return square(p.first) + square(p.second); }  // Euclidean norm squared


Point operator-(const Point &l, const Point &r) {
    return Point(l.first - r.first, l.second - r.second);
}

// Cross product calculations
ll cross(const Point &a, const Point &b) {
    return a.first * b.second - a.second * b.first;
}

ll cross(const Point &p, const Point &a, const Point &b) {
    return cross(a - p, b - p);
    
}

VectorPoints computeConvexHull(const vector<Point> &points) {
    int startIndex = 0;
    for(int i=1; i<points.size(); i++){
        if(points[startIndex].second > points[i].second || (points[startIndex].second == points[i].second && points[startIndex].first < points[i].first)){
            startIndex = i;
        }
    }

    vector<int> candidates, hull{startIndex};

    // Collecting indices of all points excluding the starting point
    for (int i = 0; i < points.size(); i++) {
        if (points[i] != points[startIndex]) {
            candidates.push_back(i);
        }
    }

    // Sort the points by polar angle wrt starting point
    sort(candidates.begin(), candidates.end(), [&](int a, int b) {
        Point vectorA = points[a] - points[startIndex];
        Point vectorB = points[b] - points[startIndex];
        ll crossProduct = cross(vectorA, vectorB);
        return crossProduct != 0 ? crossProduct > 0 : norm(vectorA) < norm(vectorB);
    });

    // Construct the convex hull
    for (const auto &candidate : candidates) {
        while (hull.size() > 1 &&
               cross(points[hull[hull.size() - 2]], points[hull.back()], points[candidate]) <= 0) {
            hull.pop_back();  // Removing points that do not form a counter-clockwise turn
        }
        hull.push_back(candidate);
    }

    VectorPoints convex_hull;
    for(auto i:hull){
        convex_hull.push_back(points[i]);
    }

    return convex_hull;
}


int main() {

    freopen("input.txt", "r", stdin);
    freopen("output.txt", "w", stdout);

        int numPoints;
        cin >> numPoints;

        if (numPoints == 0){
            cout<< 0 <<'\n';
            return 0;
        } 

        vector<Point> points(numPoints);
        for (int i = 0; i < numPoints; i++) {
            long long x, y;
            cin >> x >> y;
            points[i] = Point(x, y);
        }

        double start = omp_get_wtime();

        VectorPoints hull = computeConvexHull(points);

        double end = omp_get_wtime();

        cout<< "Final Convex Hull Points:\n";
        cout << hull.size() << '\n';
        for (const auto &index : hull) {
            cout << index.first << " " << index.second << '\n';
        }
        cout<<fixed<<"Time taken for sequential computation:- "<<end-start<<endl;
    return 0;
}
