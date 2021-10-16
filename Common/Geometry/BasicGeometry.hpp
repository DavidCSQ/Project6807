#pragma once
#include <Eigen/Dense>
#include <vector>

namespace geometry {
    template <typename T>
    using Vector3 = Eigen::Matrix<T, 3, 1>;

    // the plane is represented by (x - _p) /dot _normal = 0
    template <typename T>
    class Plane {
    public:
        Plane(Vector3<T> p, Vector3<T> normal) {
            _p = p;
            _normal = normal;
            _normal.normalize();
        }

        Vector3<T>& p() { return _p; }
        Vector3<T>& normal() { return _normal; }
        
        // return if the point is on plane
        // also fill parameter dist as the signed distance from point to plane
        bool onPlane(Vector3<T> point, T& dist) {
            dist = (point - _p).dot(_normal);
            if (std::fabs(dist) < 1e-6) {
                return true;
            } else {
                return false;
            }
        }

    private:
        Vector3<T> _p;
        Vector3<T> _normal;
    };

    template <typename T>
    class Triangle {
    public:
        Triangle(Vector3<T> v0, Vector3<T> v1, Vector3<T> v2) {
            _vertices[0] = v0;
            _vertices[1] = v1;
            _vertices[2] = v2;
        }

        Vector3<T>* vertices() { return _vertices; }
        Vector3<T>& vertices(int idx) { return _vertices[idx]; }

        
        // TODO: HW1
        // part 2.1
        // Implement the function to do intersection between triangle and plane p
        // Input: plane p
        // Output: intersections points with three edges
        // Hint:
        //      - enumerate three edges of the triangle and do intersection individually
        //      - consider the case that no intersection
        //      - consider how to avoid repeated intersection points in returned list
        bool indv(Vector3<T> v0, Vector3<T> v1,Plane<T> p, Vector3<T>& tbd){
            if (((v0.dot(p.normal())-v1.dot(p.normal()))==0)&&((p.p().dot(p.normal())-v1.dot(p.normal()))==0)){
                return false;
            }
            T lam = (p.p().dot(p.normal())-v1.dot(p.normal()))/(v0.dot(p.normal())-v1.dot(p.normal()));
            if ((lam>=0)&&(lam<=1)){
                tbd = lam*v0+(1-lam)*v1;
                return true;
            }
            return false;
        }
        bool Vertex_Cmp(Vector3<T> A, Vector3<T> B) {
            if ((A(0) < B(0) - 2e-6)||(A(0) > B(0) + 2e-6)) {
                return true;
            } else if((A(1) < B(1) - 2e-6)||(A(1) > B(1) + 2e-6)){
                return true;
            } else if((A(2) < B(2) - 2e-6)||(A(2) > B(2) + 2e-6)){
                return true;
            } else{
                return false;
            }
        }
        std::vector<Vector3<T>> IntersectPlane(Plane<T> p) {
            std::vector<Vector3<T>> intersections;
            intersections.clear();
            Vector3<T> A1;
            Vector3<T> A2;
            Vector3<T> A3;
            T dist;
            if(p.onPlane(_vertices[0],dist)&&p.onPlane(_vertices[1],dist)&&p.onPlane(_vertices[2],dist)){//if the triangle is on the plane, ignore it 
                return intersections;
            }
            bool a = indv(_vertices[0],_vertices[1],p,A1);
            bool b = indv(_vertices[1],_vertices[2],p,A2);
            bool c = indv(_vertices[2],_vertices[0],p,A3);
            if(a){//make sure that there are always exactly two points stored
                if((b)&&(Vertex_Cmp(A1,A2))){
                    intersections.push_back(A1);
                    intersections.push_back(A2);
                } else if((c)&&(Vertex_Cmp(A1,A3))){
                    intersections.push_back(A1);
                    intersections.push_back(A3);
                }
            }else if(b){
                if((c)&&(Vertex_Cmp(A2,A3))){
                    intersections.push_back(A2);
                    intersections.push_back(A3);
                }
            }
            return intersections;
        }

        // TODO: HW2
        // part 1.1
        // Implement the function to do intersection between triangle and a ray
        // Input: a ray, the ray is represented by an origin position and a direction vector
        // Output: return a real number t, the intersection is origin + dir * t, t = -1 means no intersection
        const T IntersectRay(const Vector3<T>& origin, const Vector3<T>& dir) const {
            const T flag = static_cast<T>(-1.0);    
            return flag;
        }

    private:
        Vector3<T> _vertices[3];
    };
}