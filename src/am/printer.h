#ifndef Printer_h
#define Printer_h

#include <cmath>
#include <iostream>
#include <memory>
#include <string>
#include <vector>
#include <queue>

#include "parser.h"
#include "polygon_clipping.h"

#include "floatarray.h"

// 3rdparty
#include "parallel_hashmap/phmap.h"

struct PrintMove {
    double start_position[3];
    double position[3];
    double feed_rate = 0.0;

    //double duration = 0.0;

    double extruded_length = 0.0;
    double width = 0.48;
    double height = 0.3;

    double time_start = 0.0;
    double time_end = 0.0;

    PrintMove(double x,double y,double z) {
        position[0] = x;
        position[1] = y;
        position[2] = z;
    }

    double duration() {
        return time_end - time_start;
    }

    // Calculates distance to XYZ point
    double distanceTo(double (&to)[3]) {
        return std::sqrt(
            std::pow(position[0] - to[0], 2) +
            std::pow(position[1] - to[1], 2) +
            std::pow(position[2] - to[2], 2)
        );
    }

    double dx() {
        return position[0] - start_position[0];
    }

    double dy() {
        return position[1] - start_position[1];
    }

    double dz() {
        return position[2] - start_position[2];
    }

    double length() {
        return std::sqrt(
            std::pow(position[0] - start_position[0], 2) +
            std::pow(position[1] - start_position[1], 2) +
            std::pow(position[2] - start_position[2], 2)
        );
    }

	friend std::ostream& operator<<(std::ostream& os, PrintMove& pm) {
		os << "Ts=" << pm.time_start << ", Te=" << pm.time_end << ", h=" << pm.height << ", w=" << pm.width << "L=" << pm.length();
    	return os;
	}
};

struct CellNode
{
	/// Cell index
	int index;

	// Id used for OOFEM
	int id;

	double coords[3];

	/**
	 * @brief Operator to compare two cells by index.
	 * @param o Cell
	 * @return	True	if cells are identical
	 * @return False	if cells indexes doesn't match
	*/
	bool operator==(const CellNode& o) const
	{
		return index == o.index;
	}

	/**
	 * @brief Hash Cell for usage inside of a hash map.
	 * @param p Cell to hash
	 * @return Hash
	*/
	friend size_t hash_value(const CellNode& p)
	{
		return p.index;
	}
};

struct Cell
{
	/// Cell index
	int index;

	// Id used for OOFEM
	int id;

	int nodes[8];

	/**
	 * @brief Operator to compare two cells by index.
	 * @param o Cell
	 * @return	True	if cells are identical
	 * @return False	if cells indexes doesn't match
	*/
	bool operator==(const Cell& o) const
	{
		return index == o.index;
	}

	/**
	 * @brief Hash Cell for usage inside of a hash map.
	 * @param p Cell to hash
	 * @return Hash
	*/
	friend size_t hash_value(const Cell& p)
	{
		return p.index;
	}
};

class Grid {
private:
	// Cell data (elements)
    phmap::parallel_flat_hash_map<int, Cell> map;

	// Cell nodes data (nodes)
	phmap::parallel_flat_hash_map<int, CellNode> nodes_map;

    int sizes[3] = {0, 0, 0};
	double steps[3] = {0.0, 0.0, 0.0};
	double mins[3] = {0.0, 0.0, 0.0};

public:
    Grid() {};

	/**
	 * @brief Constructs Grid with equal step size in all directions.
	 *
	 * @param step  Step size in each direction
     * @param step  Number of steps in each direction
	 */
	Grid(double step, double size)
	{
		this->steps[0] = step;
		this->steps[1] = step;
		this->steps[2] = step;

		this->sizes[0] = size;
		this->sizes[1] = size;
		this->sizes[2] = size;
	}


	/**
	 * @brief Construct Grid with given step sizes in x, y, z directions.
	 *
	 * @param x  Step to be set in x-direction
	 * @param y  Step to be set in y-direction
	 * @param z  Step to be set in z-direction
	 */
	Grid(double x, double y, double z, double nx, double ny, double nz)
	{
		this->steps[0] = x;
		this->steps[1] = y;
		this->steps[2] = z;

        this->sizes[0] = nx;
		this->sizes[1] = ny;
		this->sizes[2] = nz;
	}

	/**
	 * @brief Sets equal step size in all direction.
	 *
	 * @param step  Step to be set
	 */
	void setStep(double step)
	{
		this->steps[0] = step;
		this->steps[1] = step;
		this->steps[2] = step;
	}

	/**
	 * @brief Sets step sizes in x, y, z directions.
	 *
	 * @param x  Step to be set in x-direction
	 * @param y  Step to be set in y-direction
	 * @param z  Step to be set in z-direction
	 */
	void setStep(double x, double y, double z)
	{
		this->steps[0] = x;
		this->steps[1] = y;
		this->steps[2] = z;
	}

    void setSizes(int x, int y, int z)
	{
		this->sizes[0] = x;
		this->sizes[1] = y;
		this->sizes[2] = z;
	}

	/**
	 * @brief Returns location index given by x, y, z position.
	 *
	 * @param x Position in x-direction
	 * @param y Position in y-direction
	 * @param z Position in z-direction		      
	 * @return Calculated index
	 */
	int getIndex(int x, int y, int z)
	{
		if (x < 0 || x >= this->sizes[0])
			return -1;

		if (y < 0 || y >= this->sizes[1])
			return -1;

		if (z < 0 || z >= this->sizes[2])
			return -1;

		return x + this->sizes[0] * y + this->sizes[0] * this->sizes[1] * z;
	}

	/**
	 * @brief Returns map size.
	 *   
	 * @return Map size 
	 */
	int getSize()
	{
		return this->map.size();
	}

    /**
	 * @brief Decides whether the cell is active or not.
	 *
	 * @param i %Cell index		      
	 * @return True if cell is active     
	 * @return False if cell is inactive
	 */
	bool isActive(int i)
	{
		return !(this->map.find(i) == this->map.end());
	}

	void getVoxelGeometry(int ix, int iy, int iz, point2D (&points)[4]) {
		double x1 = ix * this->steps[0];
		double y1 = iy * this->steps[1];
		double z1 = iz * this->steps[2];

		double x2 = x1 + this->steps[0];
		double y2 = y1 + this->steps[1];
		double z2 = z1 + this->steps[2];

		points[0] = { x1, y1 };
		points[1] = { x2, y1 };
		points[2] = { x2, y2 };
		points[3] = { x1, y2 };
	}

    void add(PrintMove &p)
	{
		int newStepVoxels = 0;
        double len = p.length();
        double normal[2] = {-p.dy() / len, p.dx() / len};
        double w2 = p.width/2;

		if(std::abs(p.dx()) < 1e-12 && std::abs(p.dy()) < 1e-12)
			return;//std::cout << normal[0] << ", " << normal[1] << "\n";

        // XY polygon
        point2D subjectPolygon[] = {
            {p.start_position[0] + normal[0] * w2, p.start_position[1] + normal[1] * w2},
            {p.start_position[0] - normal[0] * w2, p.start_position[1] - normal[1] * w2},
            {p.position[0] - normal[0] * w2, p.position[1] - normal[1] * w2},
            {p.position[0] + normal[0] * w2, p.position[1] + normal[1] * w2},
        };
        int subjectPolygonSize = sizeof(subjectPolygon) / sizeof(subjectPolygon[0]);
    
		int minx = std::floor((std::min({subjectPolygon[0].x, subjectPolygon[1].x, subjectPolygon[2].x, subjectPolygon[3].x})-1e-12) / this->steps[0])-1;
		int miny = std::floor((std::min({subjectPolygon[0].y, subjectPolygon[1].y, subjectPolygon[2].y, subjectPolygon[3].y})-1e-12) / this->steps[1])-1;
		int maxx = std::ceil((std::max({subjectPolygon[0].x, subjectPolygon[1].x, subjectPolygon[2].x, subjectPolygon[3].x})+1e-12) / this->steps[0])+1;
		int maxy = std::ceil((std::max({subjectPolygon[0].y, subjectPolygon[1].y, subjectPolygon[2].y, subjectPolygon[3].y})+1e-12) / this->steps[1])+1;

        // clipping polygon
        //point2D clipPolygon[] = { {0,0}, {1000,0}, {1000,1000}, {0,1000} };
        for(int i = minx; i < maxx; i++) {
			for(int j = miny; j < maxy; j++) {
				// check if index is valid - possibly invalid due to min/maxs
				if(this->getIndex(i,j,0) == -1) continue;

				point2D clipPolygon[4];
				this->getVoxelGeometry(i, j, 0, clipPolygon);
				/*std::cout << "elem\n";
				std::cout << clipPolygon[0].x << ", " << clipPolygon[0].y << "\n";
				std::cout << clipPolygon[1].x << ", " << clipPolygon[1].y << "\n";
				std::cout << clipPolygon[2].x << ", " << clipPolygon[2].y << "\n";
				std::cout << clipPolygon[3].x << ", " << clipPolygon[3].y << "\n";*/
				int clipPolygonSize = sizeof(clipPolygon) / sizeof(clipPolygon[0]);
			
				// define the new clipped polygon (empty)
				int newPolygonSize = 0;
				int newPolygonSize2 = 0;
				point2D newPolygon[N] = { 0 };
				point2D newPolygon2[N] = { 0 };
			
				// apply clipping
				SutherlandHodgman(subjectPolygon, subjectPolygonSize, clipPolygon, clipPolygonSize, newPolygon, newPolygonSize);
				SutherlandHodgman(clipPolygon, clipPolygonSize, subjectPolygon, subjectPolygonSize, newPolygon2, newPolygonSize2);

				if(newPolygonSize2 > newPolygonSize) {
					//std::cout << "swapping\n";
					//std::cout << "nps=" << newPolygonSize << "\n";
					//std::cout << "nps2=" << newPolygonSize2 << "\n";
					std::swap(newPolygon, newPolygon2);
					//newPolygon = newPolygon2;
					newPolygonSize = newPolygonSize2;
				}

				bool completelyInside = false;
				//std::cout << "Clipped polygon points size:" << newPolygonSize << std::endl;

				// print clipped polygon points
				if(newPolygonSize >= 3 || completelyInside)
				{
					int k = (int) std::floor((p.position[2] - 1e-12)/this->steps[2]);
					int index = this->getIndex(i,j,k);

					if(index > -1) newStepVoxels++;

					if(index == -1 || this->isActive(index)) continue;

					//std::cout << "X=" << p.position[0] << "\n";
					//std::cout << "Y=" << p.position[1] << "\n";
					//std::cout << "Z=" << p.position[2] << "\n";
					//std::cout << "z=" << 0.2 * (k+1) << "\n";

					Cell c = Cell();

					// Calculate nodal indices
					// order according to brick1ht!!!!
					//
					int nodes[8] = {
						// cube top side
						index + (this->sizes[0]) * (this->sizes[1]),
						index + this->sizes[0]  + (this->sizes[0] ) * (this->sizes[1] ),
						index + 1 + this->sizes[0] + (this->sizes[0] ) * (this->sizes[1] ),
						index + 1 + (this->sizes[0] ) * (this->sizes[1] ),
						// cube bottom side
						index,
						index + this->sizes[0],
						index + 1 + this->sizes[0],
						index + 1
					};

					for(int n = 0; n < 8; n++) {
						if(!(this->nodes_map.find(nodes[n]) == this->nodes_map.end())) {
							// Node already exist, reuse id
							c.nodes[n] = this->nodes_map[nodes[n]].id;
						} else {
							// create new node id
							CellNode node = CellNode();
							node.index = nodes[n];
							node.id = this->nodes_map.size()+1;

							// compute node coords from grid to real
							node.coords[0] = (double) (node.index % this->sizes[0]) * this->steps[0];
							node.coords[1] = (double) ((node.index / this->sizes[0]) % this->sizes[1]) * this->steps[1];
							node.coords[2] = (double) std::floor(node.index / (this->sizes[0] * this->sizes[1])) * this->steps[2];

							this->nodes_map.insert({nodes[n], node});

							c.nodes[n] = node.id;
						}
					}
					
					c.id = this->map.size() + 1;
					c.index = index;

					this->map.insert({ index, c });
				}
			}
		}

		//std::cout << "total voxels=" << newStepVoxels << "\n";
        /*for(int i = 0; i < newPolygonSize; i++)
            std::cout << "(" << newPolygon[i].x << ", " << newPolygon[i].y << ")" << std::endl;
*/


		/*int indexes[3];

		for (int i = 0; i < 3; i++)
		{
			indexes[i] = (int)floor((p.getCoords()[i] - this->mins[i]) / this->steps[i]);
		}

		// Calculate Point index
		int index = indexes[0] + this->sizes[0] * indexes[1] + this->sizes[0] * this->sizes[1] * indexes[2];

		if (this->isActive(index))
		{ // Recalculate average value in cell
			//Cell &c = this->map.find(index)->second;
			Cell& c = this->map[index];
			c.x = (c.n * c.x + p.getX()) / (c.n + 1);
			c.y = (c.n * c.y + p.getY()) / (c.n + 1);
			c.z = (c.n * c.z + p.getZ()) / (c.n + 1);
			c.n = c.n + 1;
		}
		else
		{
			Cell c = Cell();
			c.index = index;
			c.x = p.getX();
			c.y = p.getY();
			c.z = p.getZ();
			c.n = 1;

			this->map.insert({ index, c });
		}*/
	}

	const phmap::parallel_flat_hash_map<int, Cell> &getCells() const
	{
		return this->map;
	}

	const phmap::parallel_flat_hash_map<int, CellNode> &getCellNodes() const
	{
		return this->nodes_map;
	}
};

struct PrintInfo {
    double total_time = 0.0;
    double distance_moved = 0.0;
    double filament_extruded_length = 0.0;
};

class Printer final {
private:
    Grid grid = Grid();

    // Print information
    PrintInfo print_info = PrintInfo();

    // Filament diameter
    double filament_diameter = 1.75;

    // Actual origin for relative positioning
    // {X, Y, Z, E}
    double origin[4] = {0.0, 0.0, 0.0, 0.0};

    // Actual printer position in the 3D space
    double position[3] = {0.0, 0.0, 0.0};

    // Actual feed rate (speed) specified by the g-code
    double feed_rate = 0.0;

    // Movement distances can be relative or absolute
    bool movement_relative_mode = false;

    // Extruder distances can be relative or absolute
    bool extruder_relative_mode = false;

    // Elapsed print time - starts from 0, incremented as printed
    double time = 0;

	int lastMoveIndex = 0;

    // Current stack of print moves used for motion planning (speed calculation)
    std::queue<PrintMove*> motionPlannerStack;

    // Length of motion planner stack
    int MOTION_PLANNER_STACK_LENGTH = 16;

    // Vector of all moves
    std::vector<PrintMove> moves;

    void calculate_motion();

public:
    Printer();

	Printer(double stepX, double stepY, double stepZ);

	~Printer() {
		std::cout << "total nodes: " << this->grid.getCellNodes().size() << "\n";
		std::cout << "total voxels: " << this->grid.getCells().size() << "\n";
	}

	const Grid &getGrid() const
	{
		return this->grid;
	}

    void begin_print() {
        // Reset origin to zero
        this->origin[0] = 0.0;
        this->origin[1] = 0.0;
        this->origin[2] = 0.0;
        this->origin[3] = 0.0;

        // Reset position to zero
        this->position[0] = 0.0;
        this->position[1] = 0.0;
        this->position[2] = 0.0;

        // reset time
        this->time = 0.0;
    };

    void parse(std::string &file);

    void print_to_time(double time);

	void set_position(gpr::block &block);

    void add_print_move(gpr::block &block);
};
#endif