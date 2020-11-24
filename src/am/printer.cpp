#define _USE_MATH_DEFINES

#include <cmath>
#include <iostream>
#include <fstream>
#include <string>
#include <unordered_map>

#include "parser.h"
#include "printer.h"

Printer::Printer() {
    //this->grid.setStep(0.25, 0.25, 0.3); // 1mm
    //this->grid.setSizes(800, 800, 1000); // 20cm = 200
}

Printer::Printer(double stepX, double stepY, double stepZ) {
    this->grid.setStep(stepX, stepY, stepZ); // 1mm
    this->grid.setSizes((int) 250 / stepX, (int) 250 / stepY, (int) 250 / stepZ); // 20cm = 200
}

void Printer::calculate_motion() {
    this->begin_print();

    double distance = 0.0;

    // Pre fill motion planner queue
    for(int i = 0; i < this->MOTION_PLANNER_STACK_LENGTH; i++) {
        this->motionPlannerStack.push(&this->moves[i]);    
    }

    // Process the queue
    for(int i = this->MOTION_PLANNER_STACK_LENGTH; i < this->moves.size() + this->MOTION_PLANNER_STACK_LENGTH; i++) {
        PrintMove *nextMove = this->motionPlannerStack.front();
        double dist = nextMove->distanceTo(this->position);

        double duration = dist / nextMove->feed_rate;
        //std::cout << this->position[0] << ", " << this->position[1] << ", " << this->position[2] << "\n";
        //std::cout << nextMove->position[0] << ", " << nextMove->position[1] << ", " << nextMove->position[2] << "\n";
        //std::cout << dist << "\n";

        this->print_info.total_time += duration;
        this->print_info.distance_moved += dist;
        this->print_info.filament_extruded_length += nextMove->extruded_length;

        nextMove->time_start = this->print_info.total_time - duration;
        nextMove->time_end = this->print_info.total_time;

        // Set move position start
        nextMove->start_position[0] = this->position[0];
        nextMove->start_position[1] = this->position[1];
        nextMove->start_position[2] = this->position[2];

        // Set current position
        this->position[0] = nextMove->position[0];
        this->position[1] = nextMove->position[1];
        this->position[2] = nextMove->position[2];

        // Pop FIFO queue and add next move pointer
        this->motionPlannerStack.pop();

        if(i < this->moves.size()) {
            this->motionPlannerStack.push(&this->moves[i]);
        }
    }
}

void Printer::set_position(gpr::block &block) {
    std::unordered_map<char, double> data;

    for(int j = 0; j < block.size(); j++) {
        gpr::chunk chunk = block.get_chunk(j);
        if(chunk.tp() == gpr::chunk_type::CHUNK_TYPE_WORD_ADDRESS && chunk.get_address().tp() == gpr::address_type::ADDRESS_TYPE_DOUBLE) {
            data[chunk.get_word()] = chunk.get_address().double_value();
        }
        else if(chunk.tp() == gpr::chunk_type::CHUNK_TYPE_WORD_ADDRESS && chunk.get_address().tp() == gpr::address_type::ADDRESS_TYPE_INTEGER) {
            data[chunk.get_word()] = (double) chunk.get_address().int_value();
        }
    }

    if(data.find('X') != data.end()) {
        this->origin[0] = data['X'];
    }

    if(data.find('Y') != data.end()) {
        this->origin[1] = data['Y'];
    }

    if(data.find('Z') != data.end()) {
        this->origin[2] = data['Z'];
    }

    if(data.find('E') != data.end()) {
        this->origin[3] = data['E'];
    }
}

void Printer::add_print_move(gpr::block &block) {
    std::unordered_map<char, double> data;

    //std::cout << block.to_string().c_str() << "\n";
    //std::cout << "Block size: " << block.size() << "\n";
            
    for(int j = 0; j < block.size(); j++) {
        gpr::chunk chunk = block.get_chunk(j);
        if(chunk.tp() == gpr::chunk_type::CHUNK_TYPE_WORD_ADDRESS && chunk.get_address().tp() == gpr::address_type::ADDRESS_TYPE_DOUBLE) {
            data[chunk.get_word()] = chunk.get_address().double_value();
        }
        else if(chunk.tp() == gpr::chunk_type::CHUNK_TYPE_WORD_ADDRESS && chunk.get_address().tp() == gpr::address_type::ADDRESS_TYPE_INTEGER) {
            data[chunk.get_word()] = (double) chunk.get_address().int_value();
        }
    }

    double feed_rate = this->feed_rate;
    double extruded_length = 0.0;

    if(data.find('X') != data.end()) {
        this->position[0] = data['X'] - this->origin[0];
    }

    if(data.find('Y') != data.end()) {
        this->position[1] = data['Y'] - this->origin[1];
    }

    if(data.find('Z') != data.end()) {
        this->position[2] = data['Z'] - this->origin[2];
    }

    if(data.find('F') != data.end()) {
        feed_rate = data['F']/60; // Set feed rate to mm/s! (F is in mm/min)
    }

    if(data.find('E') != data.end()) {
        if(this->extruder_relative_mode) {
            extruded_length = data['E'];
        } else {
            extruded_length = data['E'] - this->origin[3];
        }
    }

    PrintMove mov = PrintMove(this->position[0], this->position[1], this->position[2]);
    this->feed_rate = feed_rate;
    mov.feed_rate = feed_rate;
    mov.extruded_length = extruded_length;
    this->moves.push_back(mov);
}

void Printer::parse(std::string &file) {
    std::ifstream t(file);
    std::string file_contents((std::istreambuf_iterator<char>(t)),
			      std::istreambuf_iterator<char>());

    gpr::gcode_program p = gpr::parse_gcode(file_contents);

    std::cout << "G-code lines to parse: " << p.num_blocks() << "\n";

    // Move commands G0, G1
    gpr::chunk g0 = gpr::make_word_int('G', 0);
    gpr::chunk g1 = gpr::make_word_int('G', 1);

    // Set position
    gpr::chunk g92 = gpr::make_word_int('G', 92);

    // Absolute coordinates for extrusion
    gpr::chunk m82 = gpr::make_word_int('M', 82);

    // Relative coordinates for extrusion
    gpr::chunk m83 = gpr::make_word_int('M', 83);

    for(int i = 0; i < p.num_blocks(); i++) {
        gpr::block block = p.get_block(i);

        // Move commands
        if(block.get_chunk(0) == g1 || block.get_chunk(0) == g0) {
            this->add_print_move(block);
        }

        if(block.get_chunk(0) == g92) {
            this->set_position(block);
        }

        // Absolute extrusion
        if(block.get_chunk(0) == m82) {
            this->extruder_relative_mode = true;
        }

        // Relative extrusion
        if(block.get_chunk(0) == m83) {
            this->extruder_relative_mode = true;
        }
    }

    std::cout << "G-code moves parsed: " << this->moves.size() << "\n";

    this->calculate_motion();
    
    std::cout << "Total time: " << this->print_info.total_time << " s\n";
    std::cout << "Distance moved: " << this->print_info.distance_moved << " mm\n";
    std::cout << "Extruded filament: " << this->print_info.filament_extruded_length << " mm";
    std::cout << " | " << this->print_info.filament_extruded_length/10*M_PI*std::pow(this->filament_diameter/2/10, 2);
    std::cout << " cm^3\n";
}

void Printer::print_to_time(double target_time) {
    int nMoves = 0;
    bool finish = false;

    //this->time = 0;
    //std::cout << "Printer time " << this->time << ", target time is " << target_time << ", starting...\n";

    // loop through pre-processed moves and find those that are new
    for(int i = this->lastMoveIndex; i < this->moves.size(); i++) {
        PrintMove &move = this->moves[i];

        // Advance time
        this->time += move.duration();
        nMoves++;

        // Activate grid
        if(std::abs(move.dz()) < 1e-12 && move.extruded_length > 1e-12) {
           //std::cout << move << "\n";
           this->grid.add(move);
        }

        if(this->time > target_time) {
            this->time = std::min(target_time, move.time_end);
            this->lastMoveIndex = i + 1;
            break;
        }
    }

    std::cout << "Printer time " << this->time << ", target time was " << target_time << ", processed " << nMoves << "\n";

}