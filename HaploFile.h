
#ifndef __HAPLOFILE_H
#define __HAPLOFILE_H


#include <string>
#include <vector>

#include "Utils.h"
#include "GenoData.h"
#include "HaploPattern.h"
#include "HaploBuilder.h"


class HaploFile {
protected:
	GenoData m_genos;
	string m_filename;

	bool m_has_id;

public:
	HaploFile();
	explicit HaploFile(const string &filename);

	const string &filename() const { return m_filename; }

	void setFileName(const string &filename) { m_filename = filename; }
	void setHasID(bool enable) { m_has_id = enable; }

	virtual void readGenoData(GenoData &genos);
	virtual void writeGenoData(GenoData &genos, const char *suffix = "");

	virtual void writePattern(HaploBuilder &genos, const char *suffix = "");

	static int getFileNameNum(const string &format);
	static HaploFile *getHaploFile(const string &format, vector<string>::const_iterator fn);

protected:
	char *readAlleleName(char *buffer);
	char *writeAlleleName(char *buffer);
};

inline HaploFile::HaploFile()
: m_has_id(true)
{
}

inline HaploFile::HaploFile(const string &filename)
: m_has_id(true),
  m_filename(filename)
{
}


class HaploFileHPM : public HaploFile {
protected:
	int m_line_start;
	bool m_weighted;

public:
	HaploFileHPM() {};
	HaploFileHPM(const string &filename) : HaploFile(filename) {};

	virtual void readGenoData(GenoData &genos);
	virtual void writeGenoData(GenoData &genos, const char *suffix = NULL);
	virtual void writeGenoDataWithFreq(GenoData &genos, const char *suffix = NULL);

protected:
	void checkHeader(char *buffer);
	virtual char *readHaplotype(Haplotype &h, char *buffer);
	virtual char *writeHaplotype(const Haplotype &h, char *buffer);
	virtual void alleleTypeM2S(Allele &a);
	virtual void alleleTypeS2M(Allele &a);
};


class HaploFileHPM2 : public HaploFileHPM {
public:
	HaploFileHPM2() {};
	HaploFileHPM2(const string &filename) : HaploFileHPM(filename) {};

protected:
	virtual char *readHaplotype(Haplotype &h, char *buffer);
	virtual char *writeHaplotype(const Haplotype &h, char *buffer);
};


class HaploFileBench : public HaploFile {
protected:
	string m_children_file;
	string m_posinfo_file;

	int m_parents_num;
	int m_children_num;

public:
	HaploFileBench();
	HaploFileBench(const string &filename, const string &posinfo, const string &children = string());

	int parents_num() const { return m_parents_num; }
	int children_num() const { return m_children_num; }

	virtual void readGenoData(GenoData &genos);
	virtual void writeGenoData(GenoData &genos, const char *suffix = NULL);
	virtual void writeGenoDataWithFreq(GenoData &genos, const char *suffix = NULL);

protected:
	void readHaploFile(vector<Haplotype*> &haplos, const char *filename);
	void readPositionInfo(const char *filename);
	void writePositionInfo(const char *filename);
	char *readHaplotype(Haplotype &h, char *buffer, int heterozygous);
	char *writeHaplotype(const Haplotype &h, char *buffer);
};

inline HaploFileBench::HaploFileBench()
: m_parents_num(0),
  m_children_num(0)
{
}

inline HaploFileBench::HaploFileBench(const string &filename, const string &posinfo, const string &children)
: HaploFile(filename),
  m_children_file(children),
  m_posinfo_file(posinfo),
  m_parents_num(0),
  m_children_num(0)
{
}


#endif // __HAPLOFILE_H
