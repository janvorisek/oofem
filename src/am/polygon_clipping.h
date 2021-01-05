#ifndef Clipping
#define Clipping
 
struct point2D { double x, y; };
 
const int   N = 99; // clipped (new) polygon size
 
// check if a point is on the LEFT side of an edge
bool inside(point2D p, point2D p1, point2D p2);
 
// calculate intersection point
point2D intersection(point2D cp1, point2D cp2, point2D s, point2D e);
 
// Sutherland-Hodgman clipping
void SutherlandHodgman(point2D *subjectPolygon, int &subjectPolygonSize, point2D *clipPolygon, int &clipPolygonSize, point2D (&newPolygon)[N], int &newPolygonSize);

#endif