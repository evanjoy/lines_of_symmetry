#include <iostream>
#include <string>
#include <tuple>
#include <set>
#include <vector>
#include "./Eigen/Dense"

//using namespace std;

class SymmetrySolver {
public:
    SymmetrySolver()
    {
        rotation_by_90_ << 0, -1,
                         1, 0;
    }

    int countLinesOfSymmetry(const std::vector<Eigen::Vector2d>& points) { 
        seen_lines_.clear();
        if (points.size() < 2)
        {
            return -1;
        }
        // Assume elements in points are unique, otherwise we'd have to remove duplicates before checking size.

        int num_lines_of_symmetry = 0;
        for (int i = 0; i < points.size() - 1; i++)
        {
            for (int j = i + 1; j < points.size(); j++)
            {
                // test if the line segment (p[i], p[j]) is on a line of symmetry
                Eigen::Vector2d direction_vec = points[j] - points[i];
                if (pointsHaveLineOfSymmetry(points, direction_vec, points[i]))
                {
                    // std::cout << "*** *** Found line of symmetry along edge " << i << " to " << j << "*** ***" << std::endl;
                    ++num_lines_of_symmetry;
                }

                // test if the perpendicular bisector of line segment (p[i], p[j]) is a line of symmetry
                Eigen::Vector2d perpendicular_vec = rotation_by_90_ * direction_vec;
                if (pointsHaveLineOfSymmetry(points, perpendicular_vec, (points[i] + points[j]) / 2))
                {
                    // std::cout << "*** *** Found line of symmetry on midpoint between points " << i << " and " << j << "*** ***" << std::endl;
                    ++num_lines_of_symmetry;
                }
            }
        }
        return num_lines_of_symmetry;
    }


    /**
     * @brief Get the Ortho Direction And Distance To Origin object
     * 
     * @param line_direction normalized direction of line
     * @param line_pt a point on the line
     * @return std::pair<Eigen::Vector2d, double> the normalized perpendicular vector to the line, and how far along it to the origin
     */
    inline std::pair<Eigen::Vector2d, double> getOrthoDirectionAndDistanceToOrigin(Eigen::Vector2d& line_direction, Eigen::Vector2d line_pt)
    {
        Eigen::Vector2d perp_dir =  rotation_by_90_ * line_direction;
        double distance = perp_dir.dot(line_pt);
        return std::pair<Eigen::Vector2d, double>(perp_dir, distance);
    }


    /**
     * @brief pointsHaveLineOfSymmetry: return true if the vector of points have a line of symmetry as defined by the supplied line representation
     * 
     * @param points the vector of points
     * @param line_vec the direction of the line
     * @param line_pt a reference point along the line
     * @return true when the specified line is a line of symmetry for the 
     */
    bool pointsHaveLineOfSymmetry(const std::vector<Eigen::Vector2d>& points, const Eigen::Vector2d& line_dir, const Eigen::Vector2d& line_pt)
    {
        std::vector<Eigen::Vector2d> remaining_points(points);
        // Normalize the line direction and always make the x component positive as a convention for checking for identity / repeats
        Eigen::Vector2d normalized_dir = line_dir / line_dir.norm();
        if (normalized_dir[0] < 0 ||
            (almost_equal(normalized_dir[0], 0.0) && normalized_dir[1] < 0))
        {
            normalized_dir *= -1;
        }
        Eigen::Vector2d perp_dir;
        double distance_from_origin;
        std::tie<Eigen::Vector2d, double>(perp_dir, distance_from_origin) =
            getOrthoDirectionAndDistanceToOrigin(normalized_dir, line_pt);
        if (haveSeenLine(normalized_dir, distance_from_origin))
        {
            return false;
        }
        markLineSeen(normalized_dir, distance_from_origin);
        while (remaining_points.size() > 0)
        {
            Eigen::Vector2d pt(remaining_points.back());
            remaining_points.pop_back();
            Eigen::Vector2d relative_xy = pt - line_pt;
            double perpendicular_distance = perp_dir.dot(relative_xy);
            Eigen::Vector2d reflected_pt = pt - 2 * perpendicular_distance * perp_dir;
            if (arePointsAlmostEqual(reflected_pt, pt))
            {
                // point is on the line.  We already removed it, so continue.
                continue;
            }
            // Otherwise see if another point matches
            bool found = false;
            for (int i = 0; i < remaining_points.size(); i++)
            {
                if (arePointsAlmostEqual(remaining_points[i], reflected_pt))
                {
                    // Point matches.  Take it out too, then continue checking.
                    if (i != remaining_points.size()-1)
                    {
                        std::swap(remaining_points[i], remaining_points[remaining_points.size()-1]);
                    }
                    remaining_points.pop_back();
                    found = true;
                    break;
                }
            }
            if (!found)
            {
                return false;
            }
        }
        return true;
    }

private:
    bool haveSeenLine(Eigen::Vector2d& direction, double distance_from_origin)
    {
        // would be great if set could check almost-equal for us efficiently but we have to just do an O(n) search instead.
        //return seen_lines_.find(tuple<double, double, double>(direction[0], direction[1], distance_from_origin)) != seen_lines_.end();

        for (auto line_desc : seen_lines_)
        {
            auto& [x, y, d] = line_desc;
            if (almost_equal(x, direction[0]) && almost_equal(y, direction[1]) && almost_equal(distance_from_origin, d))
                return true;
        }
        return false;
    }

    void markLineSeen(Eigen::Vector2d& direction, double distance_from_origin)
    {
        seen_lines_.insert(std::tuple<double, double, double>(direction[0], direction[1], distance_from_origin));
    }

    void printVector(Eigen::Vector2d& vec)
    {
        std::cout << "<" << vec[0] << ", " << vec[1] << ">";
    }
 
    bool almost_equal(double x, double y)
    {
        return std::fabs(x-y) <= 1e-15;
    }

    bool arePointsAlmostEqual(Eigen::Vector2d& v1, Eigen::Vector2d& v2)
    {
        // the machine epsilon has to be scaled to the magnitude of the values used
        // and multiplied by the desired precision in ULPs (units in the last place)
        return almost_equal(v1[0], v2[0]) && almost_equal(v1[1], v2[1]);
    }

    Eigen::Matrix2d rotation_by_90_;
    std::vector<Eigen::Vector2d> midpoints_;
    std::set<std::tuple<double, double, double>> seen_lines_; // expressed as the line direction and min distance from origin
};

int main(int argc, char **argv)
{
    SymmetrySolver solver;
    std::vector<Eigen::Vector2d> points;
    Eigen::Vector2d pt;

    // 2 colinear points should have 2
    // pt << 1, 1;
    // points.push_back(pt);
    // pt << 1, 2;
    // points.push_back(pt);
    // pt << 1, 1.5;
    // points.push_back(pt);

    // 3 colinear points with 1 being centered between the others should have 2
    // pt << 1, 1;
    // points.push_back(pt);
    // pt << 1, 2;
    // points.push_back(pt);
    // pt << 1, 1.5;
    // points.push_back(pt);

    // 3 colinear points with none being centered between the others should have 1
    // pt << 1, 1;
    // points.push_back(pt);
    // pt << 2, 2;
    // points.push_back(pt);
    // pt << 1.15, 1.15;
    // points.push_back(pt);

    // square should have 4
    // pt << 0, 0;
    // points.push_back(pt);
    // pt << 0, 2;
    // points.push_back(pt);
    // pt << 2, 0;
    // points.push_back(pt);
    // pt << 2, 2;
    // points.push_back(pt);

    // rectangle should have 2
    // pt << 0, 0;
    // points.push_back(pt);
    // pt << 0, 1;
    // points.push_back(pt);
    // pt << 2, 0;
    // points.push_back(pt);
    // pt << 2, 1;
    // points.push_back(pt);

    // horizontally symmetric triangle should have 1
    // pt << 1, 1;
    // points.push_back(pt);
    // pt << 0, 0;
    // points.push_back(pt);
    // pt << 2, 0;
    // points.push_back(pt);

    // equilateral triangle should have 3
    pt << 0, 0;
    points.push_back(pt);
    pt << 0.5, sqrt(3.0)/2.0;
    points.push_back(pt);
    pt << 1, 0;
    points.push_back(pt);

    std::cout << solver.countLinesOfSymmetry(points) << std::endl;
    return 0;
}
