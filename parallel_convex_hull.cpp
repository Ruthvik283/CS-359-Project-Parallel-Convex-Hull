#include <iostream>
#include <vector>
#include <algorithm>
#include <omp.h>
#include <utility>
#include <cmath>

using namespace std;

using ll = long long;
using Point = pair<ll, ll>;
using VectorPoints = vector<Point>;

//Functions for cross product calculations
ll cross(const Point &a, const Point &b) {
    return a.first * b.second - a.second * b.first;
}
ll cross(const Point &p, const Point &fr, const Point &to) {
    return cross(Point(fr.first - p.first, fr.second - p.second), 
                 Point(to.first - p.first, to.second - p.second));
}


Point operator-(const Point &l, const Point &r) {
    return Point(l.first - r.first, l.second - r.second);
}
Point operator+(const Point &l, const Point &r) {
    return Point(l.first + r.first, l.second + r.second);
}


ll square(ll a) { return a * a; }
ll norm(const Point &p) { return square(p.first) + square(p.second); }  // Euclidean norm squared
ll norm(const Point &l, const Point &r) { return (square(r.first - l.first) + square(r.second - l.second)); }// Norm square

// Function to find the lower tangent line between two convex hulls
pair<int, int> findLowerllangent(const VectorPoints &hullA, const VectorPoints &hullB) {
    int left_point = 0, right_point = 0;
    int hullA_size = hullA.size(), hullB_size = hullB.size();

    //Finding the right most point in the left convex hull
    for(int i=1; i<hullA_size; i++){
        if(hullA[i].first > hullA[left_point].first){
            left_point = i;
        }
    }

    //Finding the left most point in the right convex hull
    for(int i=1; i<hullB_size; i++){
        if(hullB[i].first < hullB[right_point].first){
            right_point = i;
        }
    }

    while(true){
        bool move_made = 0;
        
        // Checking for clockwise turn (moving right down) 
        int next_right = (right_point + 1)%hullB_size;
        ll cross_r = cross(hullA[left_point], hullB[right_point], hullB[next_right]);
        if(cross_r < 0 || (cross_r == 0 && hullB[next_right].second < hullB[right_point].second)){
            right_point = next_right;
            move_made = 1;
        }

        // Checking for anticlockwise turn (moving left down)
        int next_left = (left_point - 1 + hullA_size)%hullA_size;
        ll cross_l = cross(hullB[right_point], hullA[left_point], hullA[next_left]);
        if(cross_l > 0 || (cross_l == 0 && hullA[next_left].second < hullA[left_point].second)){
            left_point = next_left;
            move_made = 1;
        }

        // If not turn is made, we reached the lowermost tangent
        if(!move_made){
            break;
        }
    }

    return {left_point, right_point};
}

// Function to find the upper tangent line between two convex hulls
pair<int, int> findUpperllangent(const VectorPoints &hullA, const VectorPoints &hullB) {
    int left_point = 0, right_point = 0;
    int hullA_size = hullA.size(), hullB_size = hullB.size();

    //Finding the right most point in the left convex hull
    for(int i=1; i<hullA_size; i++){
        if(hullA[i].first > hullA[left_point].first){
            left_point = i;
        }
    }

    //Finding the left most point in the right convex hull
    for(int i=1; i<hullB_size; i++){
        if(hullB[i].first < hullB[right_point].first){
            right_point = i;
        }
    }

    while(true){
        bool move_made = 0;
        
        // Checking for anticlockwise turn (moving right up) 
        int next_right = (right_point - 1 + hullB_size)%hullB_size;
        
        ll cross_r = cross(hullA[left_point], hullB[right_point], hullB[next_right]);
        if(cross_r > 0 || (cross_r == 0 && hullB[next_right].second > hullB[right_point].second)){
            right_point = next_right;
            move_made = 1;
        }

        // Checking for clockwise turn (moving left up)
        int next_left = (left_point + 1)%hullA_size;

        ll cross_l = cross(hullB[right_point], hullA[left_point], hullA[next_left]);
        if(cross_l < 0 || (cross_l == 0 && hullA[next_left].second > hullA[left_point].second)){
            left_point = next_left;
            move_made = 1;
        }

        // If no turn is made, we reached the lowermost tangent
        if(!move_made){
            break;
        }
    }
    return {left_point, right_point};
}

// Function to merge two convex hulls using tangents
VectorPoints mergeConvexHulls(const VectorPoints &hullA, const VectorPoints &hullB) {
    auto lower = findLowerllangent(hullA, hullB);
    auto upper = findUpperllangent(hullA, hullB);
    int hullA_size = hullA.size(), hullB_size = hullB.size();

    VectorPoints tempHull;
    
    // Adding all the right hull points
    for(int i=lower.second; ; i++){
        i %= hullB_size;
        tempHull.push_back(hullB[i]);
        if(i == upper.second){
            break;
        }
    }

    // Adding all the left hull points
    for(int i=upper.first; ; i++){
        i %= hullA_size;
        tempHull.push_back(hullA[i]);
        if(i == lower.first){
            break;
        }
    }

    VectorPoints mergedHull;
    if(tempHull.size() > 2){
        int startIndex = 0;
        for(int i=1; i<tempHull.size(); i++){
            if(tempHull[startIndex].second > tempHull[i].second || (tempHull[startIndex].second == tempHull[i].second && tempHull[startIndex].first < tempHull[i].first)){
                startIndex = i;
            }
        }

        int sz = tempHull.size();
        for(int i=startIndex;i<startIndex+sz;i++){
            int curr=i%sz,left=(i-1+sz)%sz,right=(i+1)%sz;
            // if(norm(tempHull[left], tempHull[curr]) + norm(tempHull[curr], tempHull[right]) != norm(tempHull[left], tempHull[right])){
            //     mergedHull.push_back(tempHull[curr]);
            // }

            // Handling the collinear points case
            if((tempHull[left]==tempHull[curr]&&mergedHull.size()>1&&tempHull[left]!=mergedHull.back())||cross(tempHull[left],tempHull[curr],tempHull[right])!=0){
                mergedHull.push_back(tempHull[curr]);
            }

        }
    }
    else mergedHull = tempHull;
        
    return mergedHull;
}

// Function to find indices of points that form the convex hull
VectorPoints computeConvexHull(const vector<Point> &points) {
    int startIndex = 0;
    for(int i=1; i<points.size(); i++){
        if(points[startIndex].second > points[i].second || (points[startIndex].second == points[i].second && points[startIndex].first < points[i].first)){
            startIndex = i;
        }
    }
    vector<int> candidates, hull{startIndex};

    // Collect indices of all points except the starting point
    for (int i = 0; i < points.size(); i++) {
        if (points[i] != points[startIndex]) {
            candidates.push_back(i);
        }
    }

    // Sort the points by polar angle with respect to the starting point
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
            hull.pop_back();  // Remove points that do not form a counter-clockwise turn
        }
        hull.push_back(candidate);
    }

    VectorPoints convex_hull;
    for(auto i:hull){
        convex_hull.push_back(points[i]);
    }

    return convex_hull;
}


void mergeSortRecursive(VectorPoints& v, long long left, long long right) {
   if (left < right) {
      if (right-left >= 32) {
         long long mid = (left+right)/2; 
         #pragma omp taskgroup
         {
            #pragma omp task shared(v) untied if(right-left >= (1<<14))
            mergeSortRecursive(v, left, mid);
            #pragma omp task shared(v) untied if(right-left >= (1<<14))
            mergeSortRecursive(v, mid+1, right);
            #pragma omp taskyield
         }
         inplace_merge(v.begin()+left, v.begin()+mid+1, v.begin()+right+1);
      }else{
         sort(v.begin()+left, v.begin()+right+1);
     }
    }
}


void mergeSort(VectorPoints& v) { 
     #pragma omp parallel
     #pragma omp single
     mergeSortRecursive(v, 0, v.size()-1); 
}

int main() {
    
    freopen("input.txt", "r", stdin);
    freopen("output.txt", "w", stdout);

        int numPoints;
        cin >> numPoints;

        VectorPoints points(numPoints);
        // Enter the points (x y) format"
        for (int i = 0; i < numPoints; ++i) {
            cin >> points[i].first >> points[i].second;
        }


        //sorting points and removing duplicates from input data
        mergeSort(points);

        double start = omp_get_wtime();
        auto last = std::unique(points.begin(), points.end());
        points.erase(last, points.end());

        numPoints = points.size();

        if(numPoints == 0){
            cout << "0\n";
            return 0;
        }


        // Number of threads
        int numThreads = 4;
        numThreads = min(numPoints, numThreads);

        // Divide points into chunks for each thread
        vector<VectorPoints> localHulls(numThreads);
        int chunkSize = numPoints / numThreads;

        #pragma omp parallel for num_threads(numThreads)
        for(int i=0;i<numThreads;i++){
            int threadID = i;
            int start = threadID * chunkSize;
            int end = (threadID == numThreads - 1) ? numPoints : start + chunkSize;

            VectorPoints cur(points.begin() + start, points.begin() + end);

            // Compute convex hull for a chunk of points
            localHulls[threadID] = computeConvexHull(cur);
        }

        // Sequentially merge the computed hulls
        while (localHulls.size() > 1) {
            int sz=localHulls.size();
            vector<VectorPoints> mergedHulls((sz+1)/2);

            #pragma omp parallel for num_threads(numThreads / 2)
            for (int i = 0; i < localHulls.size(); i += 2) {
                if (i + 1 < localHulls.size()) {
                    // Merge two adjacent hulls
                    VectorPoints mergedHull = mergeConvexHulls(localHulls[i], localHulls[i + 1]);
                    mergedHulls[i/2]=mergedHull;
                } else {
                    mergedHulls[i/2]=localHulls[i];
                }
            }

            

            localHulls.swap(mergedHulls);
        }

        // Final merged convex hull
        VectorPoints finalHull = localHulls[0];

        double end = omp_get_wtime();


        cout<< "Final Convex Hull Points:\n";
        cout<< finalHull.size() << '\n';
        for (const auto &point : finalHull) {
            cout<< point.first << " " << point.second << '\n';
        }

        cout<<fixed<<"Time taken for parallel computation:- "<<end-start<<endl;

    return 0;
}