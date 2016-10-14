#include "interface.h"

interface::interface()
{

}

void interface::writeToFile(vector<Celestial*> bodies, string filename)
{
    ofstream m_file;
    m_file.open(filename, std::ios::app);
    m_file << bodies.size() << "\n";
    m_file << "some comment" << "\n";
    for(Celestial *celestial : bodies) {

        m_file << celestial->body_radius << " " << celestial->body_name << " " << celestial->r[0] << " " << celestial->r[1] << " " << celestial->r[2] << "\n";
    }
}
