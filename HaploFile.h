
#ifndef __HAPLOFILE_H
#define __HAPLOFILE_H


#include "Utils.h"
#include "HaploData.h"
#include "HaploPattern.h"
#include "HaploBuilder.h"


class HaploFile {
protected:
	HaploData *m_haplo_data;
	char m_filename[STR_LEN_FILE_NAME];

	bool m_has_id;

public:
	HaploFile();
	explicit HaploFile(const char *filename);

	const char *filename() const { return m_filename; }

	void setFileName(const char *filename) { strcpy(m_filename, filename); }
	void setHasID(bool enable) { m_has_id = enable; }

	virtual void readHaploData(HaploData &hd);
	virtual void writeHaploData(HaploData &hd, const char *suffix = NULL);

	virtual void writePattern(HaploBuilder &hd, const char *suffix = NULL);

protected:
	char *readAlleleName(char *buffer);
	char *writeAlleleName(char *buffer);

	char *makeFileName(char *filename, const char *suffix = NULL);
};


class HaploFileHPM : public HaploFile {
protected:
	int m_line_start;
	bool m_weighted;

public:
	HaploFileHPM() {};
	HaploFileHPM(const char *filename) : HaploFile(filename) {};

	virtual void readHaploData(HaploData &hd);
	virtual void writeHaploData(HaploData &hd, const char *suffix = NULL);
	virtual void writeHaploDataWithFreq(HaploData &hd, const char *suffix = NULL);

protected:
	void checkHeader(char *buffer);
	char *readHaplotype(Haplotype &h, char *buffer);
	char *writeHaplotype(const Haplotype &h, char *buffer);
};


class HaploFileBench : public HaploFile {
protected:
	char m_children_file[STR_LEN_FILE_NAME];
	char m_posinfo_file[STR_LEN_FILE_NAME];

	int m_parents_num;
	int m_children_num;

public:
	HaploFileBench();
	HaploFileBench(const char *filename, const char *posinfo, const char *children = NULL);

	int parents_num() const { return m_parents_num; }
	int children_num() const { return m_children_num; }

	virtual void readHaploData(HaploData &hd);
	virtual void writeHaploData(HaploData &hd, const char *suffix = NULL);
	virtual void writeHaploDataWithFreq(HaploData &hd, const char *suffix = NULL);

protected:
	void readHaploFile(List<Haplotype> &haplos, const char *filename);
	char *readHaplotype(Haplotype &h, char *buffer, int heterozygous);
	char *writeHaplotype(const Haplotype &h, char *buffer);
};


#endif // __HAPLOFILE_H
