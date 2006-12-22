
#include <string.h>

#include "Allele.h"


////////////////////////////////
//
// class AlleleSequence

AlleleSequence::AlleleSequence()
{
	m_alleles = NULL;
	m_length = 0;
}

AlleleSequence::AlleleSequence(int len)
{
	int i;
	if (len > 0) {
		m_length = len;
		m_alleles = new Allele [m_length];
		for (i=0; i<m_length; i++) {
			m_alleles[i] = -1;
		}
	}
	else {
		m_alleles = NULL;
		m_length = 0;
	}
}

AlleleSequence::AlleleSequence(const AlleleSequence &as)
{
	int i;
	m_length = as.m_length;
	if (m_length > 0) {
		m_alleles = new Allele [m_length];
		for (i=0; i<m_length; i++) {
			m_alleles[i] = as.m_alleles[i];
		}
	}
	else {
		m_alleles = NULL;
		m_length = 0;
	}
}

AlleleSequence::~AlleleSequence()
{
	delete[] m_alleles;
}

bool AlleleSequence::isMatch(const AlleleSequence &as, int start1, int start2, int len) const
{
	int i;
	bool match;
	if (start1+len > m_length || start2+len > as.m_length) {
		match = false;
	}
	else {
		match = true;
		for (i=0; i<len; i++) {
			if (!isMatch(as[start2+i], start1+i)) {
				match = false;
				break;
			}
		}
	}
	return match;
}

int AlleleSequence::getDiffNum(const AlleleSequence &as, int start1, int start2, int len) const
{
	int i, match;
	match = 0;
	if (start1+len <= m_length && start2+len <= as.m_length) {
		for (i=0; i<len; i++) {
			if (!isMatch(as[start2+i], start1+i)) {
				match++;
			}
		}
	}
	return match;
}

int AlleleSequence::setLength(int i)
{
	if (m_length != i) {
		delete[] m_alleles;
		m_length = i;
		if (m_length > 0) {
			m_alleles = new Allele [m_length];
		}
		else {
			m_alleles = NULL;
			m_length = 0;
		}
	}
	return m_length;
}

char *AlleleSequence::readAllele(const char type, char *buffer, int &allele)
{
	char *buf;
	buf = buffer;
	buf += strspn(buf, " \t\r\n");					// skip delimiters
	if (type == 'S') {							// SNP, bi-allelic
		if (buf[0] == '-' || buf[0] == '?') {	// missing allele
			allele = -1;
		}
		else {
			allele = (unsigned char) buf[0];
		}
		buf++;
	}
	else {										// microsatellite
		if (buf[0] == '-' || buf[0] == '?') {	// missing allele
			allele = -1;
		}
		else {
			allele = atoi(buf);
		}
		buf += strcspn(buf, " \t\r\n");
	}
	return buf;
}

char *AlleleSequence::writeAllele(const char type, char *buffer, int &allele)
{
	char buf[64];
	buffer[0] = 0;
	if (type == 'S') {							// SNP, bi-allelic
		if (allele < 0) {						// missing allele
			buffer[0] = '?';
		}
		else {
			buffer[0] = (char) allele;
		}
		buffer[1] = ' ';							// put delimiter
		buffer[2] = 0;
	}
	else {										// microsatellite
		sprintf(buf, "%d ", allele);
		strcat(buffer, buf);
	}
	return buffer;
}

char *AlleleSequence::read(const char *types, char *buffer, int len)
{
	int i;
	char *buf;
	Allele a;
	buf = buffer;
	if (types == NULL) {
		setLength(0);
		while (buf[0] != 0 && (len <=0 || length() < len)) {
			buf = readAllele('M', buf, a);
			(*this) += a;
		};
	}
	else {
		i = strlen(types);
		if (len > 0) {
			if (i < len) {
				Logger::error("Incomplete allele type information!");
				exit(1);
			}
		}
		else {
			len = i;
		}
		setLength(len);
		for (i=0; i<m_length; i++) {
			buf = readAllele(types[i], buf, m_alleles[i]);
		}
	}
	return buffer;
}

char *AlleleSequence::write(const char *types, char *buffer) const
{
	int i;
	char *buf;
	buffer[0] = 0;
	buf = buffer;
	if (types == NULL) {
		for (i=0; i<m_length; i++) {
			writeAllele('M', buf, m_alleles[i]);
			buf = buffer + strlen(buffer);
		}
	}
	else {
		for (i=0; i<m_length; i++) {
			writeAllele(types[i], buf, m_alleles[i]);
			buf = buffer + strlen(buffer);
		}
	}
	return buffer;
}

AlleleSequence &AlleleSequence::assign(const AlleleSequence &as)
{
	int i;
	if (m_length != as.m_length) {
		delete[] m_alleles;
		m_length = as.m_length;
		m_alleles = new Allele [m_length];
	}
	for (i=0; i<m_length; i++) {
		m_alleles[i] = as.m_alleles[i];
	}
	return *this;
}

AlleleSequence &AlleleSequence::concatenate(const AlleleSequence &as)
{
	int i;
	Allele *temp = m_alleles;
	m_alleles = new Allele [m_length+as.m_length];
	for (i=0; i<m_length; i++) {
		m_alleles[i] = temp[i];
	}
	for (i=0; i<as.m_length; i++) {
		m_alleles[m_length+i] = as.m_alleles[i];
	}
	m_length += as.m_length;
	delete[] temp;
	return *this;
}

AlleleSequence &AlleleSequence::concatenate(const Allele &a)
{
	int i;
	Allele *temp = m_alleles;
	m_alleles = new Allele [m_length+1];
	for (i=0; i<m_length; i++) {
		m_alleles[i] = temp[i];
	}
	m_alleles[m_length] = a;
	m_length++;
	delete[] temp;
	return *this;
}

AlleleSequence &AlleleSequence::concatenate(const AlleleSequence &as1, const AlleleSequence &as2)
{
	int i;
	if (m_length != as1.m_length+as2.m_length) {
		delete[] m_alleles;
		m_length = as1.m_length+as2.m_length;
		m_alleles = new Allele [m_length];
	}
	for (i=0; i<as1.m_length; i++) {
		m_alleles[i] = as1.m_alleles[i];
	}
	for (i=0; i<as2.m_length; i++) {
		m_alleles[as1.m_length+i] = as2.m_alleles[i];
	}
	return *this;
}

AlleleSequence &AlleleSequence::concatenate(const AlleleSequence &as, const Allele &a)
{
	int i;
	if (m_length != as.m_length+1) {
		delete[] m_alleles;
		m_length = as.m_length+1;
		m_alleles = new Allele [m_length];
	}
	for (i=0; i<as.m_length; i++) {
		m_alleles[i] = as.m_alleles[i];
	}
	m_alleles[as.m_length] = a;
	return *this;
}

AlleleSequence &AlleleSequence::concatenate(const Allele &a, const AlleleSequence &as)
{
	int i;
	if (m_length != as.m_length+1) {
		delete[] m_alleles;
		m_length = as.m_length+1;
		m_alleles = new Allele [m_length];
	}
	m_alleles[0] = a;
	for (i=1; i<m_length; i++) {
		m_alleles[i] = as.m_alleles[i-1];
	}
	return *this;
}
