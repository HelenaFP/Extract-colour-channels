/*----------------------------------------------------------------------------------------

	Image Processing
	Assignment 3
	Department of Bioengineering
	Imperial College London
	
	

----------------------------------------------------------------------------------------*/

#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <vector>
#include <math.h>
#include "puff.h"

using namespace std;


/*----------------------------------------------------------------------------------------

	crc * adler
	
----------------------------------------------------------------------------------------*/


/* Table of CRCs of all 8-bit messages. */
unsigned long crc_table[256];

/* Flag: has the table been computed? Initially false. */
int crc_table_computed = 0;

/* Make the table for a fast CRC. */
void make_crc_table(void)
{
	unsigned long c;
 	int n, k;

 	for (n = 0; n < 256; n++)
 	{
   		c = (unsigned long) n;
   		for (k = 0; k < 8; k++)
   		{
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
  crc() routine below). */

unsigned long update_crc(unsigned long crc, unsigned char *buf, int len)
{
	unsigned long c = crc;
	int n;

	if (!crc_table_computed)
		make_crc_table();
	
	for (n = 0; n < len; n++)
	{
		c = crc_table[(c ^ buf[n]) & 0xff] ^ (c >> 8);
	}
	return c;
}

/* Return the CRC of the bytes buf[0..len-1]. */
unsigned long crc(unsigned char *buf, int len)
{
	return update_crc(0xffffffffL, buf, len) ^ 0xffffffffL;
}

const uint32_t MOD_ADLER = 65521;

uint32_t adler32(unsigned char *data, size_t len) 
{
    uint32_t a = 1, b = 0;
    size_t index;
    
    // Process each byte of the data in order
    for (index = 0; index < len; ++index)
    {
        a = (a + data[index]) % MOD_ADLER;
        b = (b + a) % MOD_ADLER;
    }
    
    return (b << 16) | a;
}

/*----------------------------------------------------------------------------------------

	struct Header
	
----------------------------------------------------------------------------------------*/

struct Header
{
	int w;			// width
	int h;			// height
	int b;			// bit depth
	int c;			// color type
	int cm;			// compression method
	int fm;			// filter method
	int im;			// interlacing method
};

/*----------------------------------------------------------------------------------------

	class Image
	
----------------------------------------------------------------------------------------*/

class Image
{
private:
	vector<vector<vector<int>>> img;
	Header head;
	
public:
	Image();
	Image(string file_name);
	~Image();
	
	void open(string file_name);
	void save_as(string file_name, int rgb_option);
	
	int get_width() {return head.w;};
	int get_height() {return head.h;};
	
	int operator()(int i, int j, int k);
	void operator= (Image& image_in);
	
private:
	void createIHDRchunk(unsigned char* carray, int* len);
	void createIDATchunk(unsigned char* carray, int* len, int rgb_option);
	void createIENDchunk(unsigned char* carray, int* len);

	void store(unsigned char* data, unsigned long len);
	void unfilter_line(unsigned char* curr_line, unsigned char* prev_line,
					unsigned char* true_line, int len, int type);
	bool cha_cmp(char* cha1, char* cha2, int len);
	unsigned int cha2uint(char* cha, int size);
	unsigned int cha2uint(unsigned char* cha, int size);
	int cha2int(char* cha, int size);
	void uint2cha(unsigned int num, char* cha, int size);
	void uint2ucha(unsigned int num, unsigned char* cha, int size);

};


/* ----- constructors and destructors ----- */

Image::Image()
{
	// do nothing - 3 dimensional vector will be empty
}

Image::Image(string file_name)
{
	open(file_name);
}

Image::~Image()
{

}

/* ----- member functions - public ----- */

void Image::open(string file_name)
{

	ifstream fin;

	char sig_ch[8];
	char len_ch[4];
	char ckt_ch[4];
	char crc_ch[4];
	
	int len = 0;
	unsigned long crc_true; 
	unsigned long crc_read;

	fin.open(file_name, ios_base::in | ios_base::binary);
	
	
	cout << endl;
	cout << "Loading " << file_name << endl;
	cout << endl;
	
	// read PNG signature
	fin.read(sig_ch, 8);
	//for (int i = 0; i < 8; i++)
	//	cout << (int)(unsigned char)sig_ch[i] << " ";
	//cout << endl;
	
	// read chunk by chunk
	while(!fin.eof())
	{
		// read a chunk
		fin.read(len_ch, 4);				// read length
		fin.read(ckt_ch, 4);				// read chunk type
		len = cha2uint(len_ch, 4);
		char* dat_ch = new char[len];
		fin.read(dat_ch, len);				// read chunk data
		fin.read(crc_ch, 4);				// read crc

		// check for data corruption
		crc_read = cha2uint(crc_ch, 4);
		unsigned char* crc_check_ch = new unsigned char[len+4];
		for (int i = 0; i < 4; i++)
			crc_check_ch[i] = ckt_ch[i];
		for (int i = 0; i < len; i++)
			crc_check_ch[i+4] = dat_ch[i];
		if (crc_read != crc(crc_check_ch, len+4))
			cout << "Data is corrupted" << endl;
		
		// display chunk type
		cout << ckt_ch << " " << len << " " << crc_read << endl;
		// cout << "crc_read: " << crc_read << endl;
		// cout << "crc_calc: " << crc(crc_check_ch, len+4) << endl;
		
		// read chunk data - specific to chunk type
		if (cha_cmp(ckt_ch, (char*)"IHDR", 4))
		{
			// extract data
			head.w =  cha2uint(&dat_ch[0], 4);
			head.h =  cha2uint(&dat_ch[4], 4);
			head.b =  (int)dat_ch[8];
			head.c =  (int)dat_ch[9];
			head.cm = (int)dat_ch[10];
			head.fm = (int)dat_ch[11];
			head.im = (int)dat_ch[12];
			
			// display data
			cout << "width:       " << head.w << endl;
			cout << "height:      " << head.h << endl;
			cout << "bitdepth:    " << head.b << endl;
			cout << "colortype:   " << head.c << endl;
			cout << "comp method: " << head.cm << endl;
			cout << "filt method: " << head.fm << endl;
			cout << "intl method: " << head.im << endl;
		}
		else if (cha_cmp(ckt_ch, (char*)"pHYs", 4))
		{
			// extract data - if needed
			// display data - if needed
		}
		else if (cha_cmp(ckt_ch, (char*)"cHRM", 4))
		{
			// extract data - if needed
			// display data - if needed
		}
		else if (cha_cmp(ckt_ch, (char*)"IDAT", 4))
		{
			int ret;
			
			unsigned char* des;
			unsigned char* src;
			unsigned long des_len;
			unsigned long src_len;
			
			// setup source data - compressed data
			// (6 bytes ommitted; flags: 2 bytes, check values: 4 bytes)
			src_len = len-6;			
			src = (unsigned char*)dat_ch+2;
			
			// setup destination data - uncompressed data
			// Each row has 1 flag bytes, W pixels, and 3 bytes per pixel (for RGB)
			// In this implementation, no alpha transparencies are needed.
			des_len = (3*(head.w)+1)*head.h;
			des = new unsigned char[des_len];
			
			//cout << "des len:     " << des_len << " bytes" << endl;
			//cout << "src len:     " << src_len << " bytes" << endl;
			
			ret = puff(des, &des_len, src, &src_len);
			
			if (ret)
				cout << "puff() failed with return code " << ret << endl;
			else
			{
				cout << "puff() succeeded uncompressing " << des_len << " bytes." << endl;
				if (src_len < des_len)
					cout << "compressed bytes unused" << endl;
			}
			
			//cout << "des len:     " << des_len << " bytes" << endl;
			//cout << "src len:     " << src_len << " bytes" << endl;
			
			store(des, des_len);
			
		}
		else if (cha_cmp(ckt_ch, (char*)"IEND", 4))
		{
			// extract data - if needed
			// display data - if needed
		}
		else
		{
			// process other headers if you'd like to implement a complete PNG reader
		}
		
		cout << endl;
		
		delete [] dat_ch;
	}
	
	fin.close();
}

void Image::save_as(string file_name, int rgb_option)
{

	ofstream fout;
	
	fout.open(file_name, ios_base::out | ios_base::binary);
	
	unsigned char  sig_ch[8] = {137,  80, 78, 71, 13, 10, 26, 10};

	unsigned char  len_ch[4];
	unsigned char* chk_ch;
	unsigned char  crc_ch[4];

	int dat_len = 0;
	unsigned long crc_val; 
	
	int len = 0;
	unsigned long crc_true; 
	unsigned long crc_read;
	
	// write PNG signature
	fout.write((char*)sig_ch,8);
	//cout << sig_ch << endl;
	//for (int i = 0; i < 8; i++)
	//	cout << (int)(unsigned char)(char)sig_ch[i] << " ";
	//cout << endl;

	// write IHDR chunk
	dat_len = 17;
	uint2ucha(13, len_ch, 4);					// create length char
	chk_ch = new unsigned char[dat_len];		// create chunk type and data char
	createIHDRchunk(chk_ch, &dat_len);
	crc_val = crc(chk_ch, 17); 					// create crc char
	uint2ucha(crc_val, crc_ch, 4);
	
	fout.write((char*)len_ch, 4);
	fout.write((char*)chk_ch, dat_len);
	fout.write((char*)crc_ch, 4);
	
	//cout << crc_val << endl;
		
	delete [] chk_ch;
	
	// write IDAT chunk
	int buffer_size = 65535;
	int num_it = ceil((3*head.w+1)*head.h/(double)buffer_size);	// # of 65,535 byte batches sent
	dat_len = (3*head.w+1)*head.h+6+5*num_it+4;	// length of chunk type and chunk name
	chk_ch = new unsigned char[dat_len];		// create chunk type and data char
	createIDATchunk(chk_ch, &dat_len, rgb_option);
	uint2ucha(dat_len-4, len_ch, 4);
	crc_val = crc(chk_ch, dat_len);			// create crc char
	uint2ucha(crc_val, crc_ch, 4);
	
	
 	// for (int i=dat_len-1000; i < dat_len; i++)
 	// 	cout << (int)chk_ch[i] << " ";
 	// cout << endl;
	
	fout.write((char*)len_ch, 4);
	fout.write((char*)chk_ch, dat_len);
	fout.write((char*)crc_ch, 4);
	
	//cout << "crc value: " << crc_val << endl;
	//cout << "crc read: " << cha2uint(crc_ch, 4) << endl;
	//cout << "num_it : " << num_it << endl;	
	
	delete [] chk_ch;
	
	// write IEND chunk
	dat_len = 4;
	uint2ucha(0, len_ch, 4);
	chk_ch = new unsigned char[dat_len];
	createIENDchunk(chk_ch, &dat_len);
	crc_val = crc(chk_ch, dat_len); 					// create crc char
	uint2ucha(crc_val, crc_ch, 4);
	
	fout.write((char*)len_ch, 4);
	fout.write((char*)chk_ch, dat_len);
	fout.write((char*)crc_ch, 4);
	
	delete [] chk_ch;
	
	fout.close();
}

void Image::createIHDRchunk(unsigned char* carray, int* len)
{
	carray[0] = (unsigned char)'I';
	carray[1] = (unsigned char)'H';
	carray[2] = (unsigned char)'D';
	carray[3] = (unsigned char)'R';
	uint2ucha(head.w, &carray[4], 4);
	uint2ucha(head.h, &carray[8], 4);
	carray[12] = (unsigned char)head.b;
	carray[13] = (unsigned char)head.c;
	carray[14] = (unsigned char)head.cm;
	carray[15] = (unsigned char)head.fm;
	carray[16] = (unsigned char)head.im;
		
	*len = 17;
}

void Image::createIDATchunk(unsigned char* carray, int* len, int rgb_option)
{

	// arrange img data into a 1-dimensional array with flags = 0
	unsigned char* img1d;
	int img1d_len = (3*head.w+1)*head.h;
	img1d = new unsigned char[img1d_len];
	int m = 0;
	for (int i=0; i<head.h; i++)
	{
		img1d[m] = 0;
		m++;
		for (int j=0; j<head.w; j++)
		{	
			for (int k=0; k<3; k++)
			{
				img1d[m] = (unsigned char)img[i][j][k];
				if (rgb_option == 1 && (k == 1 || k == 2))
				{
					img1d[m] = 0;	
				}
				else if (rgb_option == 2 && (k == 0 || k == 2))
				{
					img1d[m] = 0;	
				}
				else if (rgb_option == 3 && (k == 0 || k == 1))
				{
					img1d[m] = 0;	
				}
				
				
				m++;
			}
		}
	}
	
	// setup IDAT data (n bytes) & Adler checksum (4 bytes)
	unsigned char* data_uc;
	int buffer_size = 65535;
	int num_it = ceil(img1d_len/(double)buffer_size);
	int data_len = img1d_len+6+5*num_it;
	data_uc = new unsigned char[data_len];
	data_uc[0] = 0x78;
	data_uc[1] = 0x01;
	int n = 0;
	for (int i=0; i<num_it; i++)
	{
		if (i < (num_it-1))
		{
			data_uc[n+2]=0;
			data_uc[n+3]= (unsigned int)buffer_size>>8;
			data_uc[n+4]= (unsigned int)buffer_size;
			data_uc[n+5]= ~(unsigned char)((unsigned int)buffer_size>>8);
			data_uc[n+6]= ~(unsigned char)((unsigned int)buffer_size);
			n=n+5;
			
		}
		else
		{
			data_uc[n+2]=1;
			data_uc[n+3]= (unsigned int)(data_len-6-n+5) >> 8;
			data_uc[n+4]= (unsigned int)(data_len-6-n+5);
			data_uc[n+5]= ~((unsigned int)(data_len-6-n+5) >> 8);
			data_uc[n+6]= ~((unsigned int)(data_len-6-n+5));
			n=n+5;
		}
		
		for (int j=0; j<buffer_size; j++)
		{
			data_uc[n+2]=img1d[i*buffer_size+j];
			n++;
			
			if (n > data_len-6)
				break;
			
		}
	}
	unsigned int adl_num = adler32(img1d, img1d_len);
	uint2ucha(adl_num, &data_uc[data_len-4], 4);
	
	// include header file
	carray[0] = (unsigned char)'I';			// chunk type
	carray[1] = (unsigned char)'D';
	carray[2] = (unsigned char)'A';
	carray[3] = (unsigned char)'T';
	for (int i = 0; i < data_len; i++)		// chunk data
	{
		carray[i+4] = data_uc[i];
	}
	
	// deallocate memory
	delete [] img1d;
	delete [] data_uc;
}

void Image::createIENDchunk(unsigned char* carray, int* len)
{
	carray[0] = (unsigned char)'I';
	carray[1] = (unsigned char)'E';
	carray[2] = (unsigned char)'N';
	carray[3] = (unsigned char)'D';
	
	*len = 4;
}

int Image::operator() (int i, int j, int k)
{
	return img[i][j][k];
}

void Image::operator= (Image& image_in)
{
	img.clear();
	for (int i = 0; i < image_in.get_height(); i++)
	{
		vector<vector<int>> img_line;
		for (int j = 0; j < image_in.get_width(); j++)
		{
			vector<int> rgb;
			for (int k = 0; k < 3; k++)
			{
				rgb.push_back(image_in(i,j,k));
			}
			img_line.push_back(rgb);
		}
		img.push_back(img_line);
	}
	head = image_in.head;
}




/* ----- member functions - private ----- */

void Image::store(unsigned char* data, unsigned long len)
{
	unsigned char true_line[3*head.w+1];	// unfiltered data line;
	unsigned char curr_line[3*head.w+1];	// filtered current line - x, a
	unsigned char prev_line[3*head.w+1];	// filtered previous line - b, c
	
	int ftype = 0;
	
	for (int i = 0; i < head.h; i++)
	{
		// setup lines
		for (int j = 0; j < 3*head.w+1; j++)
		{
			if (i > 0)
				prev_line[j] = true_line[j];
			else if (j > 0)
				prev_line[j] = 0;
			true_line[j] = (unsigned int)data[i*(3*head.w+1) + j];
			curr_line[j] = (unsigned int)data[i*(3*head.w+1) + j];
		}
		
		// unfilter line
		ftype = (int)true_line[0];
		unfilter_line(curr_line, prev_line, true_line, 3*head.w+1, ftype);
		//cout << "f" << ftype << " ";

		// store line data into image data
		vector<vector<int>> img_line;
		for (int j = 1; j < 3*head.w+1; j=j+3)
		{
			vector<int> rgb;
			rgb.push_back(true_line[j]);
			rgb.push_back(true_line[j+1]);
			rgb.push_back(true_line[j+2]);
			
			img_line.push_back(rgb);
		}
		img.push_back(img_line);
	}
	
	// debugging code - to check array values
	cout << "Corner 5x5:  ";
	for (int i = 0; i < 5; i++)
	{
		if (i>0)
			cout << "             ";
		for (int j = 0; j < 5; j++)
			cout << img[i][j][0] << " " << img[i][j][1] << " " << img[i][j][2] << "\t";
		cout << endl;
	}
}

void Image::unfilter_line(unsigned char* curr_line, unsigned char* prev_line,
					unsigned char* true_line, int len, int type)
{
	if (type == 0)
	{
		// Recon(x) = Filt(x) - do nothing
	}
	else if (type ==1)
	{
		// Recon(x) = Filt(x) + Recon(a)
		unsigned int x, a;
		for (int i = 1; i < len; i++)
		{
			if (i<=3)
			{
				x = curr_line[i];
				a = 0;
			}
			else
			{
				x = curr_line[i];
				a = true_line[i-3];
			}
			true_line[i] = (x+a)%256;
		}	
	}
	else if (type ==2)
	{
		// Recon(x) = Filt(x) + Recon(b)
		unsigned int x, b;
		for (int i = 1; i < len; i++)
		{
			if (i <= 3)
			{
				x = curr_line[i];
				b = prev_line[i];
			}
			else
			{
				x = curr_line[i];
				b = prev_line[i];
			}
			
			true_line[i] = (x+b)%256;
		}	
	}
	else if (type ==3)
	{
		// Recon(x) = Filt(x) + floor((Recon(a) + Recon(b)) / 2)
		unsigned int x, a, b;
		for (int i = 1; i < len; i++)
		{
			if (i <= 3)
			{
				x = curr_line[i];
				a = 0;
				b = prev_line[i];
			}
			else
			{
				x = curr_line[i];
				a = true_line[i-3];
				b = prev_line[i];
			}
			
			true_line[i] = (unsigned int)(x + floor((int)(a+b)/(double)2))%256;
		}	
	}
	else if (type ==4)
	{
		// Recon(x) = Filt(x) + PaethPredictor(Recon(a), Recon(b), Recon(c))
		unsigned int x, a, b, c;
		unsigned int p, pa, pb, pc;
		unsigned int Pr=0;
		
		for (int i = 1; i < len; i++)
		{
			if (i <= 3)
			{
				x = curr_line[i];
				a = 0;
				b = prev_line[i];
				c = 0;
			}
			else
			{
				x = curr_line[i];
				a = true_line[i-3];
				b = prev_line[i];
				c = prev_line[i-3];
			}

			p = a+b-(int)c;
			pa = abs((int)p-(int)a);
			pb = abs((int)p-(int)b);
			pc = abs((int)p-(int)c);
			
			if ((pa <= pb) && (pa <= pc))
				Pr = a;
			else if (pb <= pc)
				Pr = b;
			else
				Pr = c;
			
			true_line[i] = (x + Pr)%256;
		}
	}	
}

bool Image::cha_cmp(char* cha1, char* cha2, int len)
{
	bool states = true;
	for (int i = 0; i < len; i++)
	{
		if (cha1[i] != cha2[i])
		{
			states = false;
			break;
		}
	}
	return states;
}

unsigned int Image::cha2uint(char* cha, int size)
{
	unsigned int num = 0;
	for (int i=0; i<size; i++)
		num = num | (unsigned char)cha[i] << 8*(4-i-1);
	return num;
}

unsigned int Image::cha2uint(unsigned char* cha, int size)
{
	unsigned int num = 0;
	for (int i=0; i<size; i++)
		num = num | (unsigned char)cha[i] << 8*(size-i-1);
	return num;
}

int Image::cha2int(char* cha, int size)
{
	int num = 0;
	for (int i=0; i<size; i++)
		num = num | (unsigned char)cha[i] << 8*(size-i-1);
	return num;
}

void Image::uint2cha(unsigned int num, char* cha, int size)
{
	for (int i=0; i<size; i++)
		cha[i] = num >> 8*(size-i-1);
}

void Image::uint2ucha(unsigned int num, unsigned char* cha, int size)
{
	for (int i=0; i<size; i++)
		cha[i] = num >> 8*(size-i-1);
}

/*----------------------------------------------------------------------------------------

	int main()
	
----------------------------------------------------------------------------------------*/

int main()
{
	string file_name;

	Image img_raw;
	Image img_r;
	
	cout << "Image Processing Software" << endl;
	cout << endl;
	cout << "Specify the name of a PNG file that you would like to process." << endl;
	cout << ">";
	
	cin >> file_name;
	
	img_raw.open(file_name+".png");

	img_raw.save_as(file_name+"_copy.png", 0);
	img_raw.save_as(file_name+"_r.png", 1);
	img_raw.save_as(file_name+"_g.png", 2);
	img_raw.save_as(file_name+"_b.png", 3);
	
	//img_r = img_raw;


	
//	Image img_r;
// 	Image img_g;
// 	Image img_b;
// 	
// 	img_r.save();
// 	img_g.save();X
// 	img_b.save();X
	

}