
#include <vector>
#include <set>

#include "HaploFile.h"

#include "MemLeak.h"


#define BUFFER_LENGTH 10000


int HaploFile::getFileNameNum(const string &format)
{
	int num = 0;
	if (format == "PHASE" || format == "HPM" || format == "HPM2") {
		num = 1;
	}
	else if (format == "BENCH2") {
		num = 2;
	}
	else if (format == "BENCH3") {
		num = 3;
	}
	return num;
}

HaploFile *HaploFile::getHaploFile(const string &format, vector<string>::const_iterator fn)
{
	HaploFile *file = 0;
	if (format == "PHASE") {
		file = new HaploFile(*fn);
	}
	else if (format == "HPM") {
		file = new HaploFileHPM(*fn);
	}
	else if (format == "HPM2") {
		file = new HaploFileHPM2(*fn);
	}
	else if (format == "BENCH2") {
		file = new HaploFileBench(*fn, *(fn+1));
	}
	else if (format == "BENCH3") {
		file = new HaploFileBench(*fn, *(fn+1), *(fn+2));
	}
	return file;
}


////////////////////////////////
//
// class HaploFile

void HaploFile::readGenoData(GenoData &genos)
{
	FILE *fp;
	char line[BUFFER_LENGTH], buf[BUFFER_LENGTH];
	char *s, *delim = " \t\r\n";
	Haplotype h1, h2;
	int i, j;
	fp = fopen(m_filename.c_str(), "r");
	if (fp == NULL) {
		Logger::error("Can not open file %s!", m_filename.c_str());
		exit(1);
	}
	fscanf(fp, "%d\n", &i);
	fscanf(fp, "%d\n", &j);
	if (i <= 0 || j <= 0) {
		Logger::error("Invalid file type!");
		exit(1);
	}
	m_genos.setGenotypeNum(i);
	m_genos.setGenotypeLen(j);
	// set loci positions
	fgets(line, BUFFER_LENGTH, fp);
	s = line + strspn(line, delim);
	if (s[0] == 'P') {
		s += strcspn(s, delim);
		s += strspn(s, delim);
		for (i=0; i<m_genos.genotype_len(); i++) {
			m_genos.setAllelePosition(i, atoi(s));
			s += strcspn(s, delim);
			s += strspn(s, delim);
		}
		fgets(line, BUFFER_LENGTH, fp);
		s = line + strspn(line, delim);
	}
	// set loci type
	for (i=0; i<m_genos.genotype_len(); i++) {
		m_genos.setAlleleType(i, s[0]);
		m_genos.setAlleleName(i, "M" + int2str(i+1));
		s++;
		s += strspn(s, delim);
	}
	for (i=0; i<m_genos.genotype_num(); i++) {
		if (m_has_id) {
			fgets(line, BUFFER_LENGTH, fp);
			if (sscanf(line, "%s", buf) > 0) {
				m_genos[i].setID(buf);
			}
		}
		else {
			m_genos[i].setID(int2str(i+1));
		}
		fgets(line, BUFFER_LENGTH, fp);
		h1.read(m_genos.allele_type().c_str(), line, m_genos.genotype_len());
		fgets(line, BUFFER_LENGTH, fp);
		h2.read(m_genos.allele_type().c_str(), line, m_genos.genotype_len());
		if (h1.length() != m_genos.genotype_len() || h2.length() != m_genos.genotype_len()) {
			Logger::error("Incorrect haplotype data for individual %d!", i);
			exit(1);
		}
		h1.setID(m_genos[i].id());
		h2.setID(m_genos[i].id());
		m_genos[i].setHaplotypes(h1, h2);
	}
	fclose(fp);
	m_genos.checkAlleleSymbol();
	genos = m_genos;
}

void HaploFile::writeGenoData(GenoData &genos, const char *suffix)
{
	FILE *fp;
	char buf[BUFFER_LENGTH], id;
	int i, j;
	m_genos = genos;
	string output_file = m_filename + suffix;
	fp = fopen(output_file.c_str(), "w");
	if (fp == NULL) {
		Logger::error("Can not open file %s!", m_filename.c_str());
		exit(1);
	}
	fprintf(fp, "%d\n", m_genos.genotype_num());
	fprintf(fp, "%d\n", m_genos.genotype_len());
	fprintf(fp, "P");
	for (i=0; i<m_genos.genotype_len(); i++) {
		fprintf(fp, " %d", m_genos.allele_postition(i));
	}
	fprintf(fp, "\n");
	fprintf(fp, "%s\n", m_genos.allele_type().c_str());
	for (i=0; i<m_genos.genotype_num(); i++) {
		id = m_genos[i].id()[0];
		if (id >= '0' && id <= '9') {
			fprintf(fp, "#%s\n", m_genos[i].id().c_str());
		}
		else {
			fprintf(fp, "%s\n", m_genos[i].id().c_str());
		}
		for (j=0; j<2; j++) {
			fprintf(fp, "%s\n", m_genos[i](j).write(m_genos.allele_type().c_str(), buf));
		}
	}
	fclose(fp);
}

void HaploFile::writePattern(HaploBuilder &hb, const char *suffix)
{
	char buf[BUFFER_LENGTH];
	m_genos = *hb.genos();
	string output_file = m_filename + suffix;
	FILE *fp = fopen(output_file.c_str(), "w");
	if (fp == NULL) {
		Logger::error("Can not open file %s!", filename);
		exit(1);
	}
	fprintf(fp, "Frequency\tLength\t%s\n", writeAlleleName(buf));
	for (int i=0; i<hb.pattern_num(); ++i) {
		const HaploPattern *hp = hb.patterns(i);
		fprintf(fp, "%f\t%d\t%s\n", hp->frequency() / m_genos.genotype_num(), hp->length(), hp->write(buf, true));
	}
	fclose(fp);
}

char *HaploFile::readAlleleName(char *buffer)
{
	int i;
	char *s, *buf, *delim = " \t\r\n";
	buf = new char [strlen(buffer)+100];
	strcpy(buf, buffer);
	s = strtok(buf, delim);
	i = 0;
	while (s != NULL && i < m_genos.genotype_len()) {
		m_genos.setAlleleName(i, s);
		s = strtok(NULL, delim);
	}
	delete[] buf;
	return buffer;
}

char *HaploFile::writeAlleleName(char *buffer)
{
	int i;
	buffer[0] = 0;
	for (i=0; i<m_genos.genotype_len(); i++) {
		strcat(buffer, m_genos.allele_name(i).c_str());
		strcat(buffer, " ");
	}
	return buffer;
}


////////////////////////////////
//
// class HaploFileHPM

void HaploFileHPM::readGenoData(GenoData &genos)
{
	FILE *fp;
	char line[BUFFER_LENGTH];
	vector<Haplotype*> haplos;
	Haplotype *h;
	int i;
	fp = fopen(m_filename.c_str(), "r");
	if (fp == NULL) {
		Logger::error("Can not open file %s!", m_filename.c_str());
		exit(1);
	}
	fgets(line, BUFFER_LENGTH, fp);
	checkHeader(line);
	while(fgets(line, BUFFER_LENGTH, fp) != NULL) {
		h = new Haplotype;
		readHaplotype(*h, line);
		if (h->length() != m_genos.genotype_len()) {
			Logger::error("Incorrect haplotype data in line %d!", haplos.size()+2);
			exit(1);
		}
		haplos.push_back(h);
	}
	if (haplos.size() % 2 != 0) {
		Logger::error("Incorrect haplotype data in line %d!", haplos.size()+2);
		exit(1);
	}
	m_genos.setGenotypeNum(haplos.size()/2);
	for (i=0; i<m_genos.genotype_num(); i++) {
		m_genos[i].setHaplotypes(*haplos[2*i], *haplos[2*i+1]);
		m_genos[i].setID(haplos[2*i]->id());
	}
	fclose(fp);
	m_genos.checkAlleleSymbol();
	genos = m_genos;
}

void HaploFileHPM::writeGenoData(GenoData &genos, const char *suffix)
{
	FILE *fp;
	char buf[BUFFER_LENGTH];
	int i, j;
	m_genos = genos;
	string output_file = m_filename + suffix;
	fp = fopen(output_file.c_str(), "w");
	if (fp == NULL) {
		Logger::error("Can not open file %s!", m_filename.c_str());
		exit(1);
	}
	fprintf(fp, "Id\t%s", writeAlleleName(buf));
	fprintf(fp, "\n");
	for (i=0; i<m_genos.genotype_num(); i++) {
		for (j=0; j<2; j++) {
			fprintf(fp, "%s\n", writeHaplotype(m_genos[i](j), buf));
		}
	}
	fclose(fp);
}

void HaploFileHPM::writeGenoDataWithFreq(GenoData &genos, const char *suffix)
{
	FILE *fp;
	char buf[BUFFER_LENGTH];
	int i, j;
	m_genos = genos;
	string output_file = m_filename + suffix;
	fp = fopen(output_file.c_str(), "w");
	if (fp == NULL) {
		Logger::error("Can not open file %s!", m_filename.c_str());
		exit(1);
	}
	fprintf(fp, "Id\t%s\tCONFIDENCE\n", writeAlleleName(buf));
	for (i=0; i<m_genos.genotype_num(); i++) {
		for (j=0; j<2; j++) {
			fprintf(fp, "%s\t%e\n", writeHaplotype(m_genos[i](j), buf), m_genos[i](j).weight());
		}
	}
	fclose(fp);
}

char *HaploFileHPM::readHaplotype(Haplotype &h, char *buffer)
{
	string allele_type(m_genos.genotype_len(), 'M');
	int i;
	char *s, *buf, *delim = " \t\r\n";
	buf = new char [strlen(buffer)+100];
	strcpy(buf, buffer);
	s = strtok(buf, delim);
	h.setID(s);
	for (i=1; i<m_line_start; i++) {
		s = strtok(NULL, delim);
	}
	s += strlen(s)+1;
	s = h.read(allele_type.c_str(), s, m_genos.genotype_len());
	if (m_weighted) {
		s += strspn(s, delim);
		h.setWeight(atof(s));
		s += strcspn(s, delim);
	}
	else {
		h.setWeight(1);
	}
	return s;
}

char *HaploFileHPM::writeHaplotype(const Haplotype &h, char *buffer)
{
	string allele_type = m_genos.allele_type();
	char *s = buffer;
	strcpy(s, h.id().c_str());
	strcat(s, "\t");
	s += strlen(s);
	h.write(allele_type.c_str(), s);
	return buffer;
}

void HaploFileHPM::checkHeader(char *buffer)
{
	int i, len;
	set<string> fields;
	vector<string> names;
	const char *delim = " \t\r\n", *s;
	char *buf;
	fields.insert("Id");
	fields.insert("Status");
	fields.insert("CONFIG_ID");
	fields.insert("CONFIDENCE");
	buf = new char [strlen(buffer)+100];
	strcpy(buf, buffer);
	s = strtok(buf, delim);
	if (s != NULL && strcmp(s, "Id") != 0) {
		Logger::error("Not a valid HPM file!");
		exit(1);
	}
	m_line_start = 1;
	s = strtok(NULL, delim);
	while (s != NULL && fields.count(s)) {
		m_line_start++;
		s = strtok(NULL, delim);
	}
	len = 0;
	while (s != NULL && !fields.count(s)) {
		len++;
		names.push_back(s);
		s = strtok(NULL, delim);
	}
	if (s != NULL && strcmp(s, "CONFIDENCE") == 0) {
		m_weighted = true;
	}
	else {
		m_weighted = false;
	}
	m_genos.setGenotypeLen(len);
	for (i=0; i<len; i++) {
		m_genos.setAlleleName(i, names[i]);
	}
	delete[] buf;
}


char *HaploFileHPM2::readHaplotype(Haplotype &h, char *buffer)
{
	string allele_type(m_genos.genotype_len(), 'S');
	int i;
	char *s, *buf, *delim = " \t\r\n";
	buf = new char [strlen(buffer)+100];
	strcpy(buf, buffer);
	s = strtok(buf, delim);
	h.setID(s);
	for (i=1; i<m_line_start; ++i) {
		s = strtok(NULL, delim);
	}
	s += strlen(s)+1;
	s = h.read(allele_type.c_str(), s, m_genos.genotype_len());
	if (m_weighted) {
		s += strspn(s, delim);
		h.setWeight(atof(s));
		s += strcspn(s, delim);
	}
	else {
		h.setWeight(1);
	}
	for (i=0; i<h.length(); ++i) {
		if (h[i] == Allele('0')) h[i] = -1;
		else if (h[i] >= Allele('1') && h[i] <= Allele('9')) h[i] = h[i].asChar() - '0';
		else if (h[i] >= Allele('A') && h[i] <= Allele('Z')) h[i] = h[i].asChar() - 'A' + 10;
		else if (h[i] >= Allele('a') && h[i] <= Allele('z')) h[i] = h[i].asChar() - 'a' + 10;
	}
	return s;
}

char *HaploFileHPM2::writeHaplotype(const Haplotype &h, char *buffer)
{
	string allele_type = m_genos.allele_type();
	int i;
	Haplotype hh = h;
	for (i=0; i<hh.length(); ++i) {
		if (allele_type[i] == 'S') {
			if (hh[i].isMissing()) hh[i] = 0;
		}
		else {
			allele_type[i] = 'S';
			if (hh[i].isMissing()) hh[i] = 0;
			else if (hh[i] >= Allele(1) && hh[i] <= Allele(9)) hh[i] = hh[i].asChar() + '0';
			else if (hh[i] >= Allele(10) && hh[i] <= Allele(35)) hh[i] = hh[i].asChar() + 'A' - 10;
		}
	}
	char *s = buffer;
	strcpy(s, hh.id().c_str());
	strcat(s, "\t");
	s += strlen(s);
	hh.write(allele_type.c_str(), s);
	return buffer;
}


////////////////////////////////
//
// class HaploFileBench

void HaploFileBench::readGenoData(GenoData &genos)
{
	FILE *fp;
	char line[BUFFER_LENGTH];
	char *s, *delim = " \t\r\n";
	vector<Haplotype*> haplos;
	int i;
	fp = fopen(m_filename.c_str(), "r");
	if (fp == NULL) {
		Logger::error("Can not open file %s!", m_filename.c_str());
		exit(1);
	}
	// get genotype length
	fgets(line, BUFFER_LENGTH, fp);
	s = line + strspn(line, delim);
	m_genos.setGenotypeLen(strcspn(s, delim));
	// set loci type
	for (i=0; i<m_genos.genotype_len(); i++) {
		m_genos.setAlleleType(i, 'S');
		m_genos.setAlleleName(i, "M" + int2str(i+1));
	}
	fclose(fp);
	readHaploFile(haplos, m_filename.c_str());
	m_parents_num = haplos.size() / 2;
	if (!m_children_file.empty()) {
		readHaploFile(haplos, m_children_file.c_str());
		m_children_num = haplos.size() / 2 - m_parents_num;
	}
	m_genos.setGenotypeNum(m_parents_num + m_children_num);
	m_genos.setUnphasedNum(m_parents_num);
	for (i=0; i<m_genos.genotype_num(); i++) {
		m_genos[i].setHaplotypes(*haplos[2*i], *haplos[2*i+1]);
		m_genos[i].setID(haplos[2*i]->id());
	}
	m_genos.checkAlleleSymbol();
	// get position info
	readPositionInfo(m_posinfo_file.c_str());
	genos = m_genos;
}

void HaploFileBench::writeGenoData(GenoData &genos, const char *suffix)
{
	FILE *fp;
	char buf[BUFFER_LENGTH];
	int i, j;
	m_genos = genos;
	string output_file = m_filename + suffix;
	fp = fopen(output_file.c_str(), "w");
	if (fp == NULL) {
		Logger::error("Can not open file %s!", m_filename.c_str());
		exit(1);
	}
	for (i=0; i<m_genos.genotype_num(); i++) {
		for (j=0; j<2; j++) {
			fprintf(fp, "%s   %d 0 %s\n", writeHaplotype(m_genos[i](j), buf), 2*i+j, m_genos[i](j).id().c_str());
		}
	}
	fclose(fp);
	writePositionInfo(m_posinfo_file.c_str());
}

void HaploFileBench::writeGenoDataWithFreq(GenoData &genos, const char *suffix)
{
	FILE *fp;
	char buf[BUFFER_LENGTH];
	int i, j;
	m_genos = genos;
	string output_file = m_filename + suffix;
	fp = fopen(output_file.c_str(), "w");
	if (fp == NULL) {
		Logger::error("Can not open file %s!", m_filename.c_str());
		exit(1);
	}
	for (i=0; i<m_genos.genotype_num(); i++) {
		for (j=0; j<2; j++) {
			fprintf(fp, "%s   %f\n", writeHaplotype(m_genos[i](j), buf), m_genos[i](j).weight());
		}
	}
	fclose(fp);
	writePositionInfo(m_posinfo_file.c_str());
}

void HaploFileBench::readHaploFile(vector<Haplotype*> &haplos, const char *filename)
{
	FILE *fp;
	char line[BUFFER_LENGTH];
	char *s, *delim = " \t\r\n";
	Haplotype *h;
	int heterozygous;
	fp = fopen(filename, "r");
	if (fp == NULL) {
		Logger::error("Can not open file %s!", filename);
		exit(1);
	}
	// read haplotypes
	heterozygous = 1;
	while(fgets(line, BUFFER_LENGTH, fp) != NULL) {
		h = new Haplotype;
		s = readHaplotype(*h, line, heterozygous);
		s += strspn(s, delim);		// next is number
		s += strcspn(s, delim);		// skip
		s += strspn(s, delim);		// next is 0
		s += strcspn(s, delim);		// skip
		s += strspn(s, delim);		// next is id
		s[strcspn(s, "\r\n")] = 0;
		h->setID(s);
		if (h->length() != m_genos.genotype_len()) {
			Logger::error("Incorrect haplotype data in line %d of %s!", haplos.size()+2, filename);
			exit(1);
		}
		haplos.push_back(h);
		heterozygous = 3 - heterozygous;
	}
	if (haplos.size() % 2 != 0) {
		Logger::error("Incorrect haplotype data in line %d!", haplos.size()+2);
		exit(1);
	}
	fclose(fp);
}

void HaploFileBench::readPositionInfo(const char *filename)
{
	FILE *fp;
	char line[BUFFER_LENGTH];
	char *s, *delim = " \t\r\n";
	int i;
	fp = fopen(filename, "r");
	if (fp == NULL) {
		Logger::error("Can not open file %s!", filename);
		exit(1);
	}
	while(fgets(line, BUFFER_LENGTH, fp) != NULL) {
		s = strtok(line, delim);
		i = atoi(s);
		if (i >= 0 && i < m_genos.genotype_len()) {
			s = strtok(NULL, delim);
			m_genos.setAlleleName(i, s);
			s = strtok(NULL, delim);
			m_genos.setAllelePosition(i, atoi(s));
		}
	}
	fclose(fp);
}

void HaploFileBench::writePositionInfo(const char *filename)
{
	FILE *fp;
	int i;
	fp = fopen(filename, "w");
	if (fp == NULL) {
		Logger::error("Can not open file %s!", filename);
		exit(1);
	}
	for (i=0; i<m_genos.genotype_len(); ++i) {
		fprintf(fp, " %d   %s   %d\n", i, m_genos.allele_name(i).c_str(), m_genos.allele_postition(i));
	}
	fclose(fp);
}

char *HaploFileBench::readHaplotype(Haplotype &h, char *buffer, int heterozygous)
{
	int i;
	char *buf, *delim = " \t\r\n";
	h.setLength(m_genos.genotype_len());
	buf = buffer + strspn(buffer, delim);
	for (i=0; i<h.length(); i++) {
		if (buf[i] == '0') {
			h[i] = -1;
		}
		else if (buf[i] == '9') {
			h[i] = '0' + heterozygous;
		}
		else {
			h[i] = (unsigned char) buf[i];
		}
	}
	buf += h.length();
	return buf;
}

char *HaploFileBench::writeHaplotype(const Haplotype &h, char *buffer)
{
	int i;

	for (i=0; i<h.length(); i++) {
		if (h[i].isMissing()) {
			buffer[i] = '0';
		}
		else {
			buffer[i] = h[i].asChar();
		}
	}
	buffer[h.length()] = 0;
	return buffer;
}
