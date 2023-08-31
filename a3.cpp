// Assignement 3: Image Processing
// Author: Helena Ferreira Pinto
// CID: 01530468
// Date: 27/11/12


#include <iostream>			// This header defines the standard input/output stream objects (for example: cin, cout).
#include <fstream>			// fstream provides templates and types to both read and write from/to files.
#include <sstream>			// sstream provides templates and types that enable operation between stream buffers and string objects.
#include <vector>			// This header defines the vector class.
#include <string>			// This header defines string types, character traits and a set of converting functions.
#include <math.h>			// This header defines a set of functions to compute common mathematical operations and transformations.
							// math.h contains the functions pow() and ceil() used in this code.

#include "puff.h"


using namespace std;		// Specifies that we will be using the standard namespace
							// which contains, for example, "cin" and "cout".


// Function prototypes
string opening_menu();
bool check_open(string file_name);
int load_image (string file_name, vector<unsigned char>* p_image);
int puff(unsigned char *dest,           /* pointer to destination pointer */
         unsigned long *destlen,        /* amount of output space */
         const unsigned char *source,   /* pointer to source data pointer */
         unsigned long *sourcelen);     /* amount of input available */		
		 
uint32_t adler32(unsigned char *data, size_t len);
void make_crc_table(void);
unsigned long update_crc(unsigned long crc, unsigned char *buf, int len);
unsigned long crc(unsigned char *buf, int len);

// Class declaration
class image_data{
private:
    // Attributes
	int IHDR_length {}, pHYs_length {}, cHRM_length {}, IEND_length {} ;	// chunk length
	unsigned int IHDR_crc {}, pHYs_crc {}, cHRM_crc {}, IEND_crc {};		// chunk crc
	int width {}, height {}, bitdepth {}, colortype {}; 					// image specifications containted in the IHDR chunk
	int comp_method {}, filt_method {}, intl_method {};
	long unsigned int IDAT_length {};
	unsigned int IDAT_crc {};
	unsigned char IHDR_chunk [25];											// full IHDR chunk, including: length, chunk type, data and crc
	unsigned char IEND_chunk[12];											// full IEND chunk, including: length, chunk type, data and crc
	unsigned long destlen {};												// length of the array that stores the uncompressed data
	unsigned char IDAT_length_output[4];									// array containing the length of IDAT chunk for output images
    vector<vector<int>> rgb_image {};										// 4 images containing, respectively, the data for rgb, red, blue and green versions of the image.
	vector<vector<int>> red_image {};										// Images are stored as vectors of vectors of ints. 
	vector<vector<int>> blue_image {};										// The innermost vector of ints represents a line.
	vector<vector<int>> green_image {};										// The outermost vector contains all the lines of the image.
	

public:
    // Methods
    vector<unsigned char> input_image {};									// input image data stored as vector of unsigned chars
    void extract_chunk_information();
	void display_image_information();
    void display_iend();
    void set_destlen();
    unsigned char* uncompression ();
    void unfilter (unsigned char* dest);
	void create_rgb_channels();
    void save_image(string filename, char channel);
    
	// Constructor
    image_data(){
    };
};

int main()
{
	string file_name {};								// name of the file to process
	vector<unsigned char> image {};
	vector<unsigned char>* p_image= &image;				// pointer to vector containing image data
	image_data brainbow;								// creates brainbow, object member of the image_data class
	unsigned char* dest {};								// pointer to array that will store the uncompressed data

	do {												
	file_name= opening_menu()+ ".png";					// formats input for filenames

    }while(!check_open(file_name));						// checks the file can be opened

	load_image (file_name, p_image);					// stores file data in vector<unsigned char> image

	for (auto i:image){									// range-based for loop that stores the image data into input_image
        brainbow.input_image.push_back(i);
	}

    brainbow.extract_chunk_information();				// extracts information from header
    brainbow.set_destlen();								// calculates the destination length based on the image dimensions
	brainbow.display_image_information();				// displays information
	
	dest= brainbow.uncompression();						// uncompresses data
	brainbow.unfilter(dest);			// unfilters data and displays corner 5x5 as a check mechanism
	brainbow.display_iend();
	
	brainbow.create_rgb_channels();						// creates raw data for red, green and blue images
	brainbow.save_image("copy.png", '0');				// creates 4 images containing, respectively:
	brainbow.save_image("brainbow_r.png", 'r');			// a copy of the original and a red, green and blue version
	brainbow.save_image("brainbow_g.png",'g');
	brainbow.save_image("brainbow_b.png", 'b');

	return 0;
}


// Opening menu
string opening_menu(){

	string file_name {""};

	cout<<"Image Processing Software\n\n";
	cout<<"Specify the name of PNG file that you would like to process." << endl;
	cout<<">";

	getline (cin, file_name, '\n');

	return file_name;
}


// Checks if the file can be opened
bool check_open(string file_name){

    ifstream input_file {};

    input_file.open(file_name);
    if(input_file.is_open())
		input_file.close();
    else {
		cout<<"\nError."<<endl;
		cout<<"We could not open ";
		cout<< file_name ;
		cout<<"."<<endl;
		cout<<"Please enter a valid file name.\n"<<endl;
		return 0;
    }

    return 1;

}

// Loads image into vector<unsigned char>
int load_image (string file_name, vector<unsigned char>* p_image){

    char c {};
    ifstream input_file {};
	input_file.open(file_name, ios_base::in|ios_base::binary);			// opens the file in binary mode
    cout<<"\nLoading "<<file_name<<endl<<endl;

	while(input_file.get(c)){
//        cout<<(unsigned char)c;
		(* p_image).push_back((unsigned char)c);
//		cout<<(unsigned int)(unsigned char)c<<endl;
	}

	input_file.close();

	return 0;
}

// Extracts information from header
// Indexes were calculated based on the chunks' structure
void image_data::extract_chunk_information(){

    int counter {};
    vector<int> crc {};													// chunk crc
    int length_data [4] {}, width_data[4] {}, height_data[4] {};

    //IHDR chunk
    for (int i=8; i<8+26;++i){											// starting at 8 because the first 8 characters are the PNG signature.
        IHDR_chunk[counter++]= input_image.at(i);						// length(4)+chunk type(4)+data(13)+crc(4)=25 (total number of chars in IHDR)
    }

    IHDR_length = (int)input_image.at(11);
    
	int cursor=16+IHDR_length;
    for (int i=cursor; i<4+cursor; ++i){
        crc.push_back((int)input_image.at(i));
    }
    for (int i=0; i<4;++i){
        IHDR_crc += ((crc.at(i))*(pow(256,(3-i))));						// converts 4 ints into a total value, effectively reading them as a byte
    }																	

    // width
    int j {};
    cursor=16;
    for (int i=cursor; i<4+cursor; ++i){
        width_data[j++]=(int)input_image.at(i);
    }
    for (int i=0; i<4;++i){
        width += ((width_data[i])*(pow(256,(3-i))));
    }

    // height
    j=0;
    cursor=20;
    for (int i=cursor; i<4+cursor; ++i){
        height_data[j++]=(int)input_image.at(i);
    }
    for (int i=0; i<4;++i){
        height += ((height_data[i])*(pow(256,(3-i))));
    }

    bitdepth=(int)input_image.at(24);
    colortype=(int)input_image.at(25);
    comp_method=(int)input_image.at(26);
    filt_method=(int)input_image.at(27);
    intl_method=(int)input_image.at(28);

	// pHYs chunk
    crc.clear();
    pHYs_length = (int)input_image.at(23+IHDR_length);
    cursor=28+IHDR_length+pHYs_length;

    for (int i=cursor; i<4+cursor; ++i){
        crc.push_back((int)input_image.at(i));
    }
    for (int i=0; i<4;++i){
        pHYs_crc += ((crc.at(i))*(pow(256,(3-i))));
    }

	// cHRM chunk
    crc.clear();
    cHRM_length=(int)input_image.at(35+IHDR_length+pHYs_length);
    cursor=40+IHDR_length+pHYs_length+cHRM_length;

    for (int i=cursor; i<4+cursor; ++i){
        crc.push_back((int)input_image.at(i));
    }
    for (int i=0; i<4;++i){
        cHRM_crc += ((crc.at(i))*(pow(256,(3-i))));
    }

	// IDAT chunk
    crc.clear();
    j=0;
    for (int i= 44+IHDR_length+pHYs_length+cHRM_length ; i< (48+IHDR_length+pHYs_length+cHRM_length);++i){
        length_data[j++]=(int)input_image.at(i);
    }
    for (int i=0; i<4;++i){
        IDAT_length += ((length_data[i])*(pow(256,(3-i))));
    }

    cursor=52+IHDR_length+pHYs_length+cHRM_length+IDAT_length;
    for (int i=cursor; i<4+ cursor; ++i){
        crc.push_back((int)input_image.at(i));
    }
    for (int i=0; i<4;++i){
        IDAT_crc += ((crc.at(i))*(pow(256,(3-i))));
    }

	// IEND chunk
    int position=0;
    cursor=56+IHDR_length+pHYs_length+cHRM_length+IDAT_length;
    for (int i=cursor; i<12+ cursor; ++i){
        IEND_chunk[position++]=input_image.at(i);
    }
    crc.clear();
    IEND_length=(int)input_image.at(59+IHDR_length+pHYs_length+IDAT_length);
    cursor=64+IHDR_length+pHYs_length+cHRM_length+IDAT_length;
    for (int i=cursor; i<4+ cursor; ++i){
        crc.push_back((int)input_image.at(i));
    }
    for (int i=0; i<4;++i){
        IEND_crc += ((crc.at(i))*(pow(256,(3-i))));
    }

}

// Uncompresses data using puff() function, which was provided
unsigned char* image_data::uncompression (){

    int puff_return {};																// return code for puff() function									
    destlen=(3*width+1)*height;														// length of the destination array, amount of output space
																					// the +1 is due to the extra byte specifying the filter type for each line
    unsigned char* dest {};           												// pointer to destination array
    unsigned long* p_destlen {};        											// pointer to length of destination array
    const unsigned char* source {};  												// pointer to source data
    unsigned long sourcelen {};														// amount of input available
    unsigned long* p_sourcelen {};     												// pointer to amount of input available


    sourcelen= IDAT_length-6;														// -6 because 6 bytes in IDAT are not raw data: 2 flags, and 4 bytes for Adler checksum
    p_sourcelen= &sourcelen;
    source= &(input_image.at(54+IHDR_length+pHYs_length+cHRM_length)); 						
    p_destlen= &destlen;
    dest= new unsigned char [destlen];												// dynamic memory allocation on the heap

    puff_return= puff(dest, p_destlen, source, p_sourcelen);

    if (!puff_return){																// puff returns 0 if it runs correctly
        cout<<"puff() succeeded uncompressing "<<destlen<<" bytes."<<endl;
    }
    else
        cout<<"Error using puff(). The return code is:"<<puff_return<<endl;

    return dest;
}

// Unfilters the data one line at a time
void image_data::unfilter(unsigned char* dest){

    int filter {};																	// filter type
    unsigned char filtered_line [3*width+1];										// initial line
    unsigned char line [3*width+1];													// unfiltered line
    unsigned char previous_line [3*width+1];										// line immediately above the current line

    for (int j=0; j<height; ++j){    												// loops through lines
        for (int i=0;i<(3*width+1);++i){  											// loops through each line, one byte at a time
            filtered_line[i]=(unsigned int)dest[j*(3*width+1)+i]; 					// updates line
            if (j!=0)
                previous_line [i]= line[i];											
            else																	// if there is no previous line, the values are assumed to be 0
                previous_line [i]= 0;
        }

        filter= (int)filtered_line[0];												// filter type is the first byte of each line
//        cout<<"filter type for line "<<j<<" is "<<filter<<endl;

        switch (filter){															// process each line according to the filter type
        case 0:{																	// no filter
        for (int i=1; i<(3*width+1);++i){
            line[i]=filtered_line[i];
        }
        }
        break;
        case 1:{																	// sub filter: Recon(x) = Filt(x) + Recon(a)
        unsigned int a {}, x {};
        for (int i=1; i<(3*width+1);++i){
            x=filtered_line[i];
            if (i<4)
                a=0;																// if there is no previous value on the line, the value is assumed to be 0
            else
                a= line[i-3];
            line[i]=(x+a)%256;														// using modulo 256 so that both the output fits into bytes
        }
        }
        break;
        case 2:{																	// up filter: Recon(x) = Filt(x) + Recon(b)
        unsigned int a {}, x {}, b{};
        for (int i=1; i<(3*width+1);++i){
            x=filtered_line[i];
            b= previous_line[i];
            line[i]=(x+b)%256;
        }
        }
        break;
        case 3:{																	// average filter: Recon(x) = Filt(x) + floor((Recon(a) + Recon(b)) / 2)
        unsigned int a {}, x {}, b{};
        for (int i=1; i<(3*width+1);++i){
            x=filtered_line[i];
            b= previous_line[i];
            if (i<4)
                a=0;
            else
                a= line[i-3];
            line[i]= ((unsigned int)(x+floor((a+b)/2)))%256;

        }
        }
        break;
        case 4:{																	// paeth filter: Recon(x) = Filt(x) + PaethPredictor(Recon(a), Recon(b), Recon(c))
        int a {}, x {}, b{};
        int c{}, pa {}, pb {}, pc {}, p {}, Pr {};
        for (int i=1; i<(3*width+1);++i){
            x=filtered_line[i];
            b= previous_line[i];
            if (i<4){
                a=0;
                c=0;
            }
            else {
                a= line[i-3];
                c= previous_line[i-3];
            }
            p = a + b -c;															// using paeth function code from the PNG specification document
            pa = abs(p - a);
            pb = abs(p - b);
            pc = abs(p - c);
            if ((pa <= pb) && (pa <= pc))
                Pr = a;
            else if (pb <= pc)
                Pr = b;
            else
                Pr = c;

            line[i]= (x+Pr)%256;
        }
        }
        break;
        }

		// storing the data
        vector<int> line_image;
        for (int i=1; i<(3*width+1); ++i){											// starting at 1 so that the filter byte is not added
            line_image.push_back((int)line[i]);										// line_image stores a line
        }
        rgb_image.push_back(line_image);											// rbg_image is a vector containing the lines
    }

    cout << "Corner 5x5:  ";														// way of checking the unfiltering is correct: display top corner
    for (int j= 0; j<5; ++j){														// lines 1 to 5
        if (j>0)
            cout << "\t\t";															// formatting output
        else
            cout<< "\t";
        for (int i = 0; i < 15; ++i){												// rbg values for pixels 1 to 5
            cout << rgb_image[j][i];
            (rgb_image[j][i]>9)? cout<<"  ": cout<<"   ";
        }
        cout << endl;
    }
}

// Creates the red, green and blue versions of the image data
void image_data::create_rgb_channels(){

    vector<int> red {}, green{}, blue {};

    for (int j=0; j<rgb_image.size(); ++j){											// loops through lines
        for (int i=0; i<(3*width+1);i+=3){											// loops through each line, one value at a time
                red.push_back(rgb_image[j][i]);
                red.push_back(0);													// green and blue values set to 0
                red.push_back(0);
        }
        red_image.push_back(red);
        red.clear();																
    }

    for (int j=0; j<rgb_image.size(); ++j){
        green.push_back(0);															// sets first red value of a line to 0
        for (int i=1; i<(3*width+1); i+=3){
                green.push_back(rgb_image[j][i]);
                green.push_back(0);
                green.push_back(0);
        }
        green_image.push_back(green);
        green.clear();
    }


    for (int j=0; j<rgb_image.size(); ++j){
        blue.push_back(0);															// sets first red and green values of a line to 0
        blue.push_back(0);
        for (int i=2; i<(3*width+1); i+=3){
                blue.push_back(rgb_image[j][i]);
                blue.push_back(0);
                blue.push_back(0);
        }
        blue_image.push_back(blue);
        blue.clear();
    }

//    cout << "\nCorner 5x5 for red image:  ";
//    for (int j= 0; j<5; ++j){
//        if (j>0)
//            cout << "\t\t\t\t";
//        else
//            cout<< "\t";
//        for (int i = 0; i < 15; ++i){
//            cout << red_image[j][i];
//            (red_image[j][i]>9)? cout<<"  ": cout<<"   ";
//        }
//        cout << endl;
//    }

//    cout << "\nCorner 5x5 for green image:  ";
//    for (int j= 0; j<5; ++j){
//        if (j>0)
//            cout << "\t\t\t\t";
//        else
//            cout<< "\t";
//        for (int i = 0; i < 15; ++i){
//            cout << green_image[j][i];
//            (green_image[j][i]>9)? cout<<"  ": cout<<"   ";
//        }
//        cout << endl;
//    }

//    cout << "\nCorner 5x5 for blue image:  ";
//    for (int j= 0; j<5; ++j){
//        if (j>0)
//            cout << "\t\t\t\t";
//        else
//            cout<< "\t";
//        for (int i = 0; i < 15; ++i){
//            cout << blue_image[j][i];
//            (blue_image[j][i]>9)? cout<<"  ": cout<<"   ";
//        }
//        cout << endl;
//    }
}

// Creates output images 
// Information about the separation of IDAT data into multiple blocks was taken from: http://calmarius.net/?lang=en&page=programming%2Fzlib_deflate_quick_reference
void image_data::save_image(string filename, char channel){

    ofstream file;
    unsigned char* output {};														// pointer to array containing IDAT data to output
    unsigned char* raw_data {};														// pointer to array containing initial IDAT (before being sliced into blocks)
    file.open(filename, ios_base:: out | ios_base:: binary);						// writing to file in binary mode
    unsigned char png_signature[8]= {137, 80, 78,71, 13, 10, 26, 10};

    // IDAT length= size image(2587607) + 6 (2 flags, 4 adler checksum) + 5*40 iterations (5 additional bits per block) = 2587813
	// converting this value to 4 bytes we get: null, 39, 124, 165
    IDAT_length_output[0]= {};
    IDAT_length_output[1]= 39;
    IDAT_length_output[2]= 124;
    IDAT_length_output[3]= 165;
	
    file.write((char*)png_signature, 8);
    file.write((char*)IHDR_chunk,25);
	file.write((char*)IDAT_length_output,4);
	
	// write IDAT data using blocks due to max size limits
    raw_data = new unsigned char[((3*width+1)*height)];							
    int k{};
    for (int j=0; j<height; ++j){													// storing the data (initially in a vector of vectors) into an array in order to write it to the file
        raw_data[k++]=0;															// adds filter type 0 for every line
        for (int i=0; i<(3*width); ++i){
            if (channel=='0')														// rgb image
                raw_data[k++]=(unsigned char) rgb_image[j][i];
            if (channel=='r')														// red image
                raw_data[k++]=(unsigned char) red_image[j][i];
            if (channel=='g')														// green image
                raw_data[k++]=(unsigned char) green_image[j][i];
            if (channel=='b')														// blue image
                raw_data[k++]=(unsigned char) blue_image[j][i];
        }
    }

    int block_size= 65535;															// maximum size of each block
	unsigned int max_byte = 255;
    int iteration = ceil(((3*width+1)*height)/(double)block_size);					// number of iterations for the total data to be outputted
    unsigned int length= ((3*width+1)*height)+ 6 + 5*iteration+8;					// length of IDAT chunk, including: chunk type, data, crc
    output= new unsigned char [length];												
    int position= 6;																// start at position 6, since 1-4 contain chunk type, and 5-6 are flag bytes
    output[0]= 'I';
    output[1]= 'D';
    output[2]= 'A';
    output[3]= 'T';
    output[4]= 0x78;																					// from specification
    output[5]= 0x01;																					// flags

    for (int i=1; i<=iteration; ++i){
        if (i<iteration){
            output[position++]=0;																		// bfinal
            output[position++]=(unsigned char)255;														// length of data block expressed in 2 bytes
            output[position++]=(unsigned char)255;														// 2 bytes of 255, thus representing 65535
            output[position++]= ~(unsigned char)max_byte;												// nlength, 2's complement of length
            output[position++]= ~(unsigned char)max_byte;
        }
        else{																							// for last iteration
            output[position++]=1;																		// bfinal set to 0 to indicate it is the last block
            output[position++]=(unsigned char)((((3*width+1)*height)-block_size*iteration)>>8);			// length= image size - 65535*number of iterations  {right shifted by 8 to divide by 256}
            output[position++]=(unsigned char)(((3*width+1)*height)-block_size*iteration);
            output[position++]=~(unsigned char)((((3*width+1)*height)-block_size*iteration)>>8);
            output[position++]=~(unsigned char)(((3*width+1)*height)-block_size*iteration);
        }

        for (int j=0; j<block_size;++j){
            output[position]= raw_data[(i-1)*block_size+j];												// adds raw data 
            ++position;
            if(position>(length-8))																		// if position> length - (4 bytes adler checksum) - (4 bytes crc)
//          if(position>2587813)
                break;
        }
    }
	
    unsigned int adler= adler32(raw_data,destlen );														// gets Adler32 checksum
	// cout<<"adler: "<<channel<<" "<<adler<<endl;
    for (int i=0; i<4; ++i){
        output[position++]=adler >> 8*(3-i);															// right shifts to convert unsigned int to 4 unsigned chars
    }

    unsigned int value_crc= crc(&output[4],(length-4));													// calculates crc value for IDAT chunk
    for (int i=0; i<4; ++i){
        output[position++]=value_crc >> 8*(4-i-1);
    }
   // cout<<position<<endl;

    file.write((char*)output,length); 																	// writes IDAT data
    file.write((char*)IEND_chunk,12);																	// writes IEND chunk

    delete [] raw_data;																					// deallocates memory
	delete [] output;

}


void image_data::set_destlen(){
    
	destlen=(3*width+1)*height;

}

void image_data::display_image_information(){

	cout<<"IHDR "<< IHDR_length<<" "<<IHDR_crc<<endl;
	cout<<"width:\t\t"<<width<<endl;
	cout<<"height:\t\t"<<height<<endl;
	cout<<"bitdepth:\t"<<bitdepth<<endl;
	cout<<"colortype:\t"<<colortype<<endl;
	cout<<"comp method:\t"<<comp_method<<endl;
	cout<<"filt method:\t"<<filt_method<<endl;
	cout<<"intl method:\t"<<intl_method<<endl<<endl;
	cout<<"pHYs "<<pHYs_length<<" "<<pHYs_crc<<endl<<endl;
	cout<<"cHRM "<<cHRM_length<<" "<<cHRM_crc<<endl<<endl;
	cout<<"IDAT "<<IDAT_length<<" "<<IDAT_crc<<endl;

}

void image_data::display_iend(){

    cout<<endl<<"IEND "<<IEND_length<<" "<<IEND_crc<<endl;

}


// with help from a GTA
// description of adler32 algorithm was taken from: https://en.wikipedia.org/wiki/Adler-32
uint32_t adler32(unsigned char *data, size_t len){	
	uint32_t modulus = 65521;							// uint32_t is a fixed width integer type
    uint32_t a = 1, b = 0;
    size_t index;

    for (index = 0; index < len; ++index)				// perform calculations one byte at a time
    {
        a = (a + data[index]) % modulus;
        b = (b + a) % modulus;
    }
    return (b << 16) | a;
}

// This block of code allows the computation of the crc.
// I didn't write it.
// It was taken from the PNG specification document: https://www.w3.org/TR/PNG/#D-CRCAppendix

/* Table of CRCs of all 8-bit messages. */
   unsigned long crc_table[256];

   /* Flag: has the table been computed? Initially false. */
   int crc_table_computed = 0;

   /* Make the table for a fast CRC. */
   void make_crc_table(void){
     unsigned long c;
     int n, k;

     for (n = 0; n < 256; n++) {
       c = (unsigned long) n;
       for (k = 0; k < 8; k++) {
         if (c & 1)
           c = 0xedb88320L ^ (c >> 1);
         else
           c = c >> 1;
       }
       crc_table[n] = c;
     }
     crc_table_computed = 1;
   }

   /* Update a running CRC with the bytes buf[0..len-1]--the CRC
      should be initialized to all 1's, and the transmitted value
      is the 1's complement of the final running CRC (see the
      crc() routine below)). */

   unsigned long update_crc(unsigned long crc, unsigned char *buf,
                            int len)
   {
     unsigned long c = crc;
     int n;

     if (!crc_table_computed)
       make_crc_table();
     for (n = 0; n < len; n++) {
       c = crc_table[(c ^ buf[n]) & 0xff] ^ (c >> 8);
     }
     return c;
   }

   /* Return the CRC of the bytes buf[0..len-1]. */
   unsigned long crc(unsigned char *buf, int len)
   {
     return update_crc(0xffffffffL, buf, len) ^ 0xffffffffL;
   }
