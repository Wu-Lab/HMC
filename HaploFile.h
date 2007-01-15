
#ifndef __HAPLOFILE_H
#define __HAPLOFILE_H


#include <string>

#include "Utils.h"
#include "HaploData.h"
#include "HaploPattern.h"
#include "HaploBuilder.h"


class HaploFile {
protected:
	HaploData *m_haplo_data;
	string m_filename;

	bool m_has_id;

public:
	HaploFile();
	explicit HaploFile(const char *filename);

	const string &filename() const { return m_filename; }

	void setFileName(const char *filename) { m_filename = filename; }
	void setHasID(bool enable) { m_has_id = enable; }

	virtual void readHaploData(HaploData &hd);
	virtual void writeHaploData(HaploData &hd, const char *suffix = NULL);

	virtual void writePattern(HaploBuilder &hd, const char *suffix = NULL);

protected:
	char *readAlleleName(char *buffer);
	char *writeAlleleName(char *buffer);
};

inline HaploFile::HaploFile()
: m_haplo_data(0),
  m_has_id(true)
{
}

inline HaploFile::HaploFile(const char *filename)
: m_haplo_data(0),
  m_has_id(true),
  m_filename(filename)
{
}


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
	string m_children_file;
	string m_posinfo_file;

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

inline HaploFileBench::HaploFileBench()
: m_parents_num(0),
  m_children_num(0)
{
}

inline HaploFileBench::HaploFileBench(const char *filename, const char *posinfo, const char *children)
: HaploFile(filename),
  m_children_file(children),
  m_posinfo_file(posinfo),
  m_parents_num(0),
  m_children_num(0)
{
}


#endif // __HAPLOFILE_H
