#ifndef __STK_IO_POINTSET__
#define __STK_IO_POINTSET__

#include <fstream>
#include <sstream>
#include <map>
#include <typeinfo>

#include <stk/pointset.hpp>
#include <stk/vector.hpp>
#include <stk/exception.hpp>
#include <stk/delaunay.hpp>

namespace stk
{

namespace io
{

class PointSetStream
{
	public:
	
		enum PositionType
		{
			POS_DOUBLE = 0,
			POS_FLOAT,
			POS_INT,
			POS_INT16,
			POS_INT8
		};
		
		enum ValueType
		{
			VAL_NONE = 0,
			VAL_DOUBLE,
			VAL_FLOAT,
			VAL_INT,
			VAL_UINT,
			VAL_UINT24,
			VAL_COMPLEXD
		};
};

/////////////////////////////////////////////////////
/////////////////////// INPUT ///////////////////////
/////////////////////////////////////////////////////

template<int DIM, typename POS, typename VAL_OLD = double>
class PointSetInputStream : public PointSetStream
{
	protected:
		enum FileType
		{
			DAT = 0,
			PTS1,
			PTS2,
			PTS3
		};
		
		std::ifstream m_file;
		Domain<DIM, POS>* m_domain;
		FileType m_fileType;
		ValueType m_valueType;
		PositionType m_positionType;
		bool m_multiPts;
		bool m_binary;
		POS m_firstElement;
		int m_counter;
	
		void readHeader()
		{
			//Dat test
			if(m_file >> m_firstElement)
			{
				m_multiPts = false;
				m_fileType = DAT;
				m_valueType = VAL_NONE;
				m_domain = new UnitDomain<DIM, POS>();
				m_domain->toroidal(true);
				m_positionType = POS_DOUBLE;
				return;
			}
			else m_file.clear();
			
			std::string version;
			m_file >> version;
			
			m_multiPts = true;
			m_fileType = PTS2;
			m_binary = false;
			
			if(version == "pts1")
			{
				m_multiPts = false;
				m_fileType = PTS1;
			}
			else if(version == "pts2")
			{
				m_multiPts = false;
			}
			else if(version == "mpts1")
			{
				m_fileType = PTS1;
			}
			else if(version == "mpts2")
			{
				
			}
			else if(version == "bpts3")
			{
				m_binary = true;
				m_fileType = PTS3;
			}
			else
			{
				throw exception::Message("Unknown pts version", STK_DBG_INFO)
					.addVar("header", version);
			}
			
			//Check dimension
			int dim = -1;
			m_file >> dim;
			if(dim != DIM)
			{
				std::stringstream sserror;
				sserror << "Wrong number of dimension in PTS file : ";
				sserror << dim << " != " << DIM;
				throw exception::Message(sserror.str(), STK_DBG_INFO);
			}
			
			//Check position type
			if(m_fileType == PTS3)
			{
				std::string pType;
				m_file >> pType;
				if(pType == "d") m_positionType = POS_DOUBLE;
				else if(pType == "f") m_positionType = POS_FLOAT;
				else if(pType == "i") m_positionType = POS_INT;
				else if(pType == "i16") m_positionType = POS_INT16;
				else if(pType == "i8") m_positionType = POS_INT8;
				else
				{
					throw exception::Message("Unknown position type", STK_DBG_INFO)
						.addVar("type", pType);
				}
			}
			else m_positionType = POS_DOUBLE;
			
			//Check weight
			std::string wType;
			m_file >> wType;
			if(wType == "n") m_valueType = VAL_NONE;
			else if(wType == "d") m_valueType = VAL_DOUBLE;
			else if(wType == "f") m_valueType = VAL_FLOAT;
			else if(wType == "c") m_valueType = VAL_COMPLEXD;
			else if(wType == "i") m_valueType = VAL_INT;
			else if(wType == "ui") m_valueType = VAL_UINT;
			else if(wType == "ui24") m_valueType = VAL_UINT24;
			else
			{
				throw exception::Message("Unknown value type", STK_DBG_INFO)
					.addVar("type", wType);
			}
			
			//Check bound
			if(m_fileType == PTS2 || m_fileType == PTS3)
			{
				std::string bType;
				m_file >> bType;
				
				if(bType == "unit" || bType == "unit-tor")
				{
					m_domain = new UnitDomain<DIM, POS>();
					m_domain->toroidal(bType == "unit-tor");
				}
				else if(bType == "rect" || bType == "rect-tor")
				{
					Vector<DIM, POS> dMin, dMax;
					for(int i=0; i<DIM; i++)
					{
						m_file >> dMin[i] >> dMax[i];
						if(dMin > dMax)
						{
							std::stringstream sserror;
							sserror << "Wrong space boundaries (min>max) : ";
							sserror << dMin << " > " << dMax;
							throw exception::Message(sserror.str(), STK_DBG_INFO);
						}
					}
				
					m_domain = new RectangularDomain<DIM, POS>(dMin, dMax);
					m_domain->toroidal(bType == "rect-tor");
				}
				else if(bType == "base" || bType == "base-tor")
				{
					std::vector< Vector<DIM, POS> > base;
					Vector<DIM, POS> ref;
					m_file >> ref;
				
					for(int i=0; i<DIM; i++)
					{
						Vector<DIM, POS> v;
						m_file >> v;
						base.push_back(v);
					}
				
					m_domain = new BaseDomain<POS>(ref, base);
					m_domain->toroidal(bType == "base-tor");
				}
				else if(bType == "hexa" || bType == "hexa-tor")
				{
					m_domain = new HexagonalDomain<POS>();
					m_domain->toroidal(bType == "hexa-tor");
				}
			}
			else
			{
				std::string bType;
				bool isToroidal;
				m_file >> bType;
				if(bType != "n" && bType != "t")
				{
					std::stringstream sserror;
					sserror << "Unknown bound type : ";
					sserror << bType;
					throw exception::Message(sserror.str(), STK_DBG_INFO);
				}
				else
				{
					isToroidal = (bType=="t");
				}
				
				//Check dim limit
				Vector<DIM, POS> dMin, dMax;
				for(int i=0; i<DIM; i++)
				{
					m_file >> dMin[i] >> dMax[i];
					if(dMin > dMax)
					{
						std::stringstream sserror;
						sserror << "Wrong space boundaries (min>max) : ";
						sserror << dMin << " > " << dMax;
						throw exception::Message(sserror.str(), STK_DBG_INFO);
					}
				}
				
				if(
					dMin == Vector<DIM, POS>(0.0) &&
					dMax == Vector<DIM, POS>(1.0)
					)
				{
					m_domain = new UnitDomain<DIM, POS>();
					m_domain->toroidal(isToroidal);
				}
				else
				{
					m_domain = new RectangularDomain<DIM, POS>(dMin, dMax);
					m_domain->toroidal(isToroidal);
				}
			}
			
			if(m_binary)
			{
				char _trash;
				m_file.read((char*)&_trash, sizeof(char));
			}
		}
		
	public:
		PointSetInputStream()
		{
			m_valueType = VAL_NONE;
			m_positionType = POS_DOUBLE;
			m_fileType = DAT;
			m_counter = 0;
			m_domain = NULL;
			m_binary = false;
		}
		
		PointSetInputStream(const std::string& i_filename)
		{
			m_valueType = VAL_NONE;
			m_positionType = POS_DOUBLE;
			m_fileType = DAT;
			m_counter = 0;
			m_domain = NULL;
			m_binary = false;
			
			open(i_filename);
		}
		
		~PointSetInputStream()
		{
			close();
			if(m_domain != NULL) delete m_domain;
		}
		
		void open(const std::string& i_filename)
		{
			m_file.open(i_filename.c_str(), std::ios::binary);
			if(!m_file) throw exception::FileNotFound(i_filename, STK_DBG_INFO);
			
			try
			{
				readHeader();
			}
			catch(stk::exception::Message& e)
			{
				throw e.addVar("filename", i_filename);
			}
		}
		
		template<typename VAL>
		void read(PointSet<DIM, POS, VAL>& o_pts)
		{
			o_pts.embedDomain(m_domain->clone());
			o_pts.clear();
			
			if(m_binary)
			{
				unsigned int nPts;
				m_file.read((char*) &nPts, sizeof(unsigned int));
				o_pts.reserve(nPts);
				
				Vector<DIM, POS> pos;
				for(unsigned int i=0; i<nPts; i++)
				{
					switch(m_positionType)
					{
						case POS_DOUBLE:
							{
								Vector<DIM, double> posRead;
								m_file.read((char*) &posRead, DIM*sizeof(double));
								for(int j=0; j<DIM; j++) pos[j] = posRead[j];
							}
							break;
						case POS_FLOAT:
							{
								Vector<DIM, float> posRead;
								m_file.read((char*) &posRead, DIM*sizeof(float));
								for(int j=0; j<DIM; j++) pos[j] = posRead[j];
							}
							break;
						case POS_INT:
							{
								Vector<DIM, int> posRead;
								m_file.read((char*) &posRead, DIM*sizeof(int));
								for(int j=0; j<DIM; j++) pos[j] = posRead[j];
							}
							break;
						case POS_INT16:
							{
								Vector<DIM, short> posRead;
								m_file.read((char*) &posRead, DIM*sizeof(short));
								for(int j=0; j<DIM; j++) pos[j] = posRead[j];
							}
							break;
						case POS_INT8:
							{
								Vector<DIM, char> posRead;
								m_file.read((char*) &posRead, DIM*sizeof(char));
								for(int j=0; j<DIM; j++) pos[j] = posRead[j];
							}
							break;
						default:
							throw exception::Message("Unknown position type", STK_DBG_INFO)
								.addVar("positionType", m_positionType);
					}
					
					switch(m_valueType)
					{
						case VAL_NONE:
							o_pts.push_back(Point<DIM, POS, VAL>(pos, o_pts.getDefaultVal()));
							break;
						case VAL_DOUBLE:
							{
								double valRead;
								m_file.read((char*) &valRead, sizeof(double));
								o_pts.push_back(Point<DIM, POS, VAL>(pos, static_cast<VAL>(valRead)));
							}
							break;
						case VAL_FLOAT:
							{
								float valRead;
								m_file.read((char*) &valRead, sizeof(float));
								o_pts.push_back(Point<DIM, POS, VAL>(pos, static_cast<VAL>(valRead)));
							}
							break;
						case VAL_INT:
							{
								int valRead;
								m_file.read((char*) &valRead, sizeof(int));
								o_pts.push_back(Point<DIM, POS, VAL>(pos, static_cast<VAL>(valRead)));
							}
							break;
						case VAL_UINT:
							{
								unsigned int valRead;
								m_file.read((char*) &valRead, sizeof(unsigned int));
								o_pts.push_back(Point<DIM, POS, VAL>(pos, static_cast<VAL>(valRead)));
							}
							break;
						case VAL_UINT24:
							{
								unsigned int valRead = 0x0;
								m_file.read((char*) &valRead, 3*sizeof(unsigned char));
								o_pts.push_back(Point<DIM, POS, VAL>(pos, static_cast<VAL>(valRead)));
							}
							break;
						case VAL_COMPLEXD:
							{
								Complexd valRead;
								m_file.read((char*) &valRead, sizeof(Complexd));
								o_pts.push_back(Point<DIM, POS, VAL>(pos, static_cast<VAL>(std::real(valRead))));
							}
							break;
						default:
							throw exception::Message("Unknown value type", STK_DBG_INFO)
								.addVar("valueType", m_valueType);
					}
				}
			}
			else if(m_fileType == DAT)
			{
				Vector<DIM, POS> pos;
				
				if(m_counter == 0)
				{
					pos[0] = m_firstElement;
					for(int i=1; i<DIM; i++) m_file >> pos[i];
					o_pts.push_back(Point<DIM, POS, VAL>(pos, o_pts.getDefaultVal()));
				}
				
				while(m_file >> pos)
				{
					o_pts.push_back(Point<DIM, POS, VAL>(pos, o_pts.getDefaultVal()));
				}
			}
			else
			{
				Vector<DIM, POS> pos;
				VAL val;
				
				while(m_file >> pos)
				{
					if(m_valueType != VAL_NONE) m_file >> val;
					else val = o_pts.getDefaultVal();
					o_pts.push_back(Point<DIM, POS, VAL>(pos, val));
				}
			}
		}
		
		void ignore()
		{
			if(m_binary)
			{
				unsigned int nPts;
				std::size_t coordSz, valSz;

				m_file.read((char*) &nPts, sizeof(unsigned int));

				switch(m_positionType)
				{
					case POS_DOUBLE:
						coordSz = sizeof(double);
						break;
					case POS_FLOAT:
						coordSz = sizeof(float);
						break;
					case POS_INT:
						coordSz = sizeof(int);
						break;
					case POS_INT16:
						coordSz = sizeof(short);
						break;
					case POS_INT8:
						coordSz = sizeof(char);
						break;
					default:
						throw exception::Message("Unknown position type", STK_DBG_INFO)
							.addVar("positionType", m_positionType);
				}

				switch(m_valueType)
				{
					case VAL_NONE:
						valSz = 0;
						break;
					case VAL_DOUBLE:
						valSz = sizeof(double);
						break;
					case VAL_FLOAT:
						valSz = sizeof(float);
						break;
					case VAL_INT:
						valSz = sizeof(int);
						break;
					case VAL_UINT:
						valSz = sizeof(unsigned int);
						break;
					case VAL_UINT24:
						valSz = 3*sizeof(unsigned char);
						break;
					case VAL_COMPLEXD:
						valSz = sizeof(Complexd);
						break;
					default:
						throw exception::Message("Unknown value type", STK_DBG_INFO)
							.addVar("valueType", m_valueType);
				}
				m_file.ignore( nPts*(DIM*coordSz + valSz) );
			}
			else
			{
				Vector<DIM, POS> pos;
				VAL_OLD val;
				while(m_file >> pos)
				{
					if(m_valueType != VAL_NONE) m_file >> val;
				}
			}

		}
		
		bool next()
		{
			if(m_binary)
			{
				m_file.peek();
				return !m_file.eof();
			}
			else
			{
				if(m_file.fail())
				{
					m_file.clear();
					std::string sep;
					m_file >> sep;
					if(sep != "#patch-separator" && sep != "")
					{
						std::stringstream sserror;
						sserror << "Wrong patch separator : " << sep << " != #patch-separator";
						throw exception::Message(sserror.str(), STK_DBG_INFO);
					}
				}
				return m_file.good();
			}
		}
		
		void close()
		{
			m_file.close();
		}
		
		ValueType valueType() const
		{
			return m_valueType;
		}
		
		PositionType positionType() const
		{
			return m_positionType;
		}
		
		Domain<DIM, POS>* domain()
		{
			return m_domain;
		}
};

template<int DIM, typename POS, typename VAL>
void read(std::string i_filename, std::vector< PointSet<DIM, POS, VAL> >& o_ptsList)
{
	stk::io::PointSetInputStream<DIM, POS, VAL> ptsStream;
	ptsStream.open(i_filename);
	
	do
	{
		o_ptsList.push_back(stk::PointSet<DIM, POS, VAL>());
		ptsStream.read(o_ptsList.back());
	}
	while(ptsStream.next());
}

template<int DIM, typename POS, typename VAL>
void read(std::string i_filename, PointSet<DIM, POS, VAL>& o_pts)
{
	stk::io::PointSetInputStream<DIM, POS, VAL> ptsStream;
	ptsStream.open(i_filename);
	ptsStream.read(o_pts);
}



//////////////////////////////////////////////////////
/////////////////////// OUTPUT ///////////////////////
//////////////////////////////////////////////////////

template<int DIM, typename POS, typename VAL>
class PointSetOutputStream : public PointSetStream
{
	protected:
		std::ofstream m_file;
		Domain<DIM, POS>* m_domain;
		bool m_multiPts;
		bool m_header;
		bool m_binary;
		PositionType m_positionType;
		ValueType m_valueType;

		bool writeHeader(const PointSet<DIM, POS, VAL>& o_pts)
		{
			if(m_domain != NULL)
			{
				if(m_binary) m_file << "bpts3" << '\t';
				else if(m_multiPts) m_file << "mpts2" << '\t';
				else m_file << "pts2" << '\t';
			}
			else m_file << "pts1" << '\t';

			m_file << DIM << '\t';
			
			//Position type
			if(m_binary)
			{
				switch(m_positionType)
				{
					case POS_DOUBLE:
						m_file << "d" << '\t';
						break;
					case POS_FLOAT:
						m_file << "f" << '\t';
						break;
					case POS_INT:
						m_file << "i" << '\t';
						break;
					case POS_INT16:
						m_file << "i16" << '\t';
						break;
					case POS_INT8:
						m_file << "i8" << '\t';
						break;
					default:
						throw exception::Message("Unknown position type", STK_DBG_INFO)
							.addVar("positionType", m_positionType);
				}
			}

			//Value Type
			switch(m_valueType)
			{
				case VAL_NONE:
					m_file << "n" << '\t';
					break;
				case VAL_DOUBLE:
					m_file << "d" << '\t';
					break;
				case VAL_FLOAT:
					m_file << "f" << '\t';
					break;
				case VAL_INT:
					m_file << "i" << '\t';
					break;
				case VAL_UINT:
					m_file << "ui" << '\t';
					break;
				case VAL_UINT24:
					m_file << "ui24" << '\t';
					break;
				case VAL_COMPLEXD:
					m_file << "c" << '\t';
					break;
				default:
					exception::Message("Unknown value type", STK_DBG_INFO)
						.addVar("valueType", m_valueType);
			}
	
			if(m_domain != NULL)
			{
				const UnitDomain<DIM, POS>* unitDomain = dynamic_cast<const UnitDomain<DIM, POS>*>(m_domain);
				const RectangularDomain<DIM, POS>* rectDomain = dynamic_cast<const RectangularDomain<DIM, POS>*>(m_domain);
				const BaseDomain<POS>* baseDomain = dynamic_cast<const BaseDomain<POS>*>(m_domain);
				const HexagonalDomain<POS>* hexaDomain = dynamic_cast<const HexagonalDomain<POS>*>(m_domain);

				if(unitDomain != NULL)
				{
					if(m_domain->isToroidal()) m_file << "unit-tor";
					else m_file << "unit";
					m_file << std::endl;
				}
				else if(rectDomain != NULL)
				{
					if(m_domain->isToroidal()) m_file << "rect-tor";
					else m_file << "rect";
					m_file << std::endl;

					for(int i=0; i<DIM; i++)
					{
						m_file << rectDomain->boundingBoxMin()[i] << '\t' << rectDomain->boundingBoxMax()[i];
						if(i<DIM-1) m_file << '\t';
						else m_file << std::endl;
					}
				}
				else if(baseDomain != NULL)
				{
					if(m_domain->isToroidal()) m_file << "base-tor";
					else m_file << "base";
					m_file << std::endl;
					
					m_file
						<< baseDomain->referencePoint()[0] << '\t'
						<< baseDomain->referencePoint()[1] << '\t'
						<< baseDomain->vector(0)[0] << '\t'
						<< baseDomain->vector(0)[1] << '\t'
						<< baseDomain->vector(1)[0] << '\t'
						<< baseDomain->vector(1)[1] << std::endl;
				}
				else if(hexaDomain != NULL)
				{
					if(m_domain->isToroidal()) m_file << "hexa-tor" << '\t';
					else m_file << "hexa";
					m_file << std::endl;
				}
				else throw exception::Message("Unknown domain", STK_DBG_INFO);
			}
			else
			{
				
				if(o_pts.isToroidal()) m_file << "t" << '\t';
				else m_file << "n" << '\t';
				m_file << std::endl;
				
				m_file << o_pts.boundingBoxMin()[0]
					<< '\t' << o_pts.boundingBoxMax()[0];
				for(int i=1; i<DIM; i++)
				{
					m_file
						<< '\t' << o_pts.boundingBoxMin()[i]
						<< '\t' << o_pts.boundingBoxMax()[i];
				}
				m_file << std::endl;
			}
			return m_file.good();
		}
		
		void writePositionBinary(const Vector<DIM, POS>& position)
		{
			switch(m_positionType)
			{
				case POS_DOUBLE:
				{
					Vector<2, double> posWrite;
					for(int d=0; d<DIM; d++) posWrite[d] = position[d];
					m_file.write((char*)&posWrite, DIM*sizeof(double));
				}
				break;
				case POS_FLOAT:
				{
					Vector<2, float> posWrite;
					for(int d=0; d<DIM; d++) posWrite[d] = position[d];
					m_file.write((char*)&posWrite, DIM*sizeof(float));
				}
				break;
				case POS_INT:
				{
					Vector<2, int> posWrite;
					for(int d=0; d<DIM; d++) posWrite[d] = position[d];
					m_file.write((char*)&posWrite, DIM*sizeof(int));
				}
				break;
				case POS_INT16:
				{
					Vector<2, short> posWrite;
					for(int d=0; d<DIM; d++) posWrite[d] = position[d];
					m_file.write((char*)&posWrite, DIM*sizeof(short));
				}
				break;
				case POS_INT8:
				{
					Vector<2, char> posWrite;
					for(int d=0; d<DIM; d++) posWrite[d] = position[d];
					m_file.write((char*)&posWrite, DIM*sizeof(char));
				}
				break;
				default:
					throw exception::Message("Unknown position type", STK_DBG_INFO)
						.addVar("positionType", m_positionType);
			}
		}
		
		void writeValueBinary(const VAL& value)
		{
			switch(m_valueType)
			{
				case VAL_NONE:
				break;
				
				case VAL_DOUBLE:
				{
					double valWrite = value;
					m_file.write((char*)&valWrite, sizeof(double));
				}
				break;
				
				case VAL_FLOAT:
				{
					float valWrite = value;
					m_file.write((char*)&valWrite, sizeof(float));
				}
				break;
				
				case VAL_INT:
				{
					int valWrite = value;
					m_file.write((char*)&valWrite, sizeof(int));
				}
				break;
				
				case VAL_UINT:
				{
					unsigned int valWrite = value;
					m_file.write((char*)&valWrite, sizeof(unsigned int));
				}
				break;
				
				case VAL_UINT24:
				{
					unsigned int valWrite = value;
					m_file.write((char*)&valWrite, 3*sizeof(unsigned char));
				}
				break;
				
				case VAL_COMPLEXD:
				{
					Complexd valWrite = Complexd(value, 0);
					m_file.write((char*)&valWrite, sizeof(Complexd));
				}
				break;
				
				default:
					throw exception::Message("Unknown value type", STK_DBG_INFO)
						.addVar("valueType", m_valueType);
			}
		}
		
		void writeSubBinary(const PointSet<DIM, POS, VAL>& i_pts)
		{
			unsigned int nPts = i_pts.size();
			m_file.write((char*)&nPts, sizeof(unsigned int));
			for(int pti=0; pti<i_pts.size(); pti++)
			{
				writePositionBinary(i_pts[pti].pos());
				writeValueBinary(i_pts[pti].val());
			}
		}
		
		void writePositionText(const Vector<DIM, POS>& position)
		{
			m_file << position;
		}
		
		void writeValueText(const VAL& value)
		{
			switch(m_valueType)
			{
				case VAL_NONE:		break;
				case VAL_DOUBLE:	m_file << static_cast<double>(value); break;
				case VAL_FLOAT:		m_file << static_cast<float>(value); break;
				case VAL_INT:		m_file << static_cast<int>(value); break;
				case VAL_UINT:		m_file << static_cast<unsigned int>(value); break;
				case VAL_UINT24:	m_file << static_cast<unsigned int>(value); break;
				case VAL_COMPLEXD:	m_file << Complexd(static_cast<double>(value), 0); break;
				default:
					throw exception::Message("Unknown value type", STK_DBG_INFO)
						.addVar("valueType", m_valueType);
			}
		}

		void writeSubText(const PointSet<DIM, POS, VAL>& i_pts)
		{
			for(int pti=0; pti<i_pts.size(); pti++)
			{
				writePositionText(i_pts.at(pti).pos());
				m_file << '\t';
				writeValueText(i_pts.at(pti).val());
				m_file << std::endl;
			}
		}
		
		bool next()
		{
			if(!m_binary)
			{
				m_file << "#patch-separator" << std::endl;
			}
			return m_file.good();
		}
		
	public:
		PointSetOutputStream(const PointSetOutputStream<DIM, POS, VAL>& a)
		{
			m_domain = NULL;
			m_multiPts = true;
			m_header = true;
			m_binary = false;
			
			if(typeid(POS) == typeid(float)) m_positionType = POS_FLOAT;
			else if(typeid(POS) == typeid(int)) m_positionType = POS_INT;
			else if(typeid(POS) == typeid(short)) m_positionType = POS_INT16;
			else if(typeid(POS) == typeid(char)) m_positionType = POS_INT8;
			else m_positionType = POS_DOUBLE;
			
			setWeight(true);
		}
		
		PointSetOutputStream()
		{
			m_domain = NULL;
			m_multiPts = true;
			m_header = true;
			m_binary = false;
			
			if(typeid(POS) == typeid(float)) m_positionType = POS_FLOAT;
			else if(typeid(POS) == typeid(int)) m_positionType = POS_INT;
			else if(typeid(POS) == typeid(short)) m_positionType = POS_INT16;
			else if(typeid(POS) == typeid(char)) m_positionType = POS_INT8;
			else m_positionType = POS_DOUBLE;
			
			setWeight(true);
		}
		
		PointSetOutputStream(const std::string& i_filename, const bool& i_multiPts = true, const bool& i_weight = true, const bool& i_header = true)
		{			
			m_domain = NULL;
			m_multiPts = true;
			m_header = true;
			m_binary = false;
			
			if(typeid(POS) == typeid(float)) m_positionType = POS_FLOAT;
			else if(typeid(POS) == typeid(int)) m_positionType = POS_INT;
			else if(typeid(POS) == typeid(short)) m_positionType = POS_INT16;
			else if(typeid(POS) == typeid(char)) m_positionType = POS_INT8;
			else m_positionType = POS_DOUBLE;
			
			setWeight(i_weight);
			
			open(i_filename, i_multiPts, i_header);
		}
		
		~PointSetOutputStream()
		{
			close();
			if(m_domain != NULL) delete m_domain;
		}

		void setMultiPts(const bool& i_multiPts = true) { m_multiPts = i_multiPts; }
		void setWeight(const bool& i_weight = true)
		{
			if(!i_weight) m_valueType = VAL_NONE;
			else
			{
				if(typeid(VAL) == typeid(double)) m_valueType = VAL_DOUBLE;
				else if(typeid(VAL) == typeid(float)) m_valueType = VAL_FLOAT;
				else if(typeid(VAL) == typeid(int)) m_valueType = VAL_INT;
				else if(typeid(VAL) == typeid(unsigned int)) m_valueType = VAL_UINT;
				else if(typeid(VAL) == typeid(Complexd)) m_valueType = VAL_COMPLEXD;
				else m_valueType = VAL_NONE;
			}
		}
		void setHeader(const bool& i_formatDat = true) { m_header = i_formatDat; }
		void setBinary(const bool& i_binary = true) { m_binary = i_binary; }
		void setValueType(const ValueType& i_valueType)
		{
			m_valueType = i_valueType;
		}
		void setPositionType(const PositionType& i_positionType)
		{
			m_positionType = i_positionType;
		}

		void open(const std::string& i_filename, const bool& i_multiPts = true, const bool& i_header = true)
		{
			m_multiPts = i_multiPts;
			m_header = i_header;
			m_file.open(i_filename.c_str(), std::ios::binary);
			if(!m_file) throw exception::FileNotWritable(i_filename, STK_DBG_INFO);
		}

		void write(const PointSet<DIM, POS, VAL>& i_pts)
		{
			if(m_domain == NULL)
			{
				m_domain = i_pts.domain()->clone();
				if(m_header) writeHeader(i_pts);
			}
			else if(m_multiPts && m_header) next();
			else throw exception::Message("Trying to write more than 1 pointset in pts/dat file, use multi-pts file (mpts) instead", STK_DBG_INFO);

			if(!m_binary) m_file.precision(15);

			if(!m_header)
			{
				for(int pti=0; pti<i_pts.size(); pti++)
				{
					m_file << i_pts.at(pti).pos() << std::endl;
				}
			}
			else if(m_binary) writeSubBinary(i_pts);
			else writeSubText(i_pts);
		}

		void close() { if(m_file.is_open()) m_file.close(); }
};


template<int DIM, typename POS, typename VAL>
void writePts(const std::string i_filename, const std::vector< PointSet<DIM, POS, VAL> >& o_ptsList, const bool& i_weight = true)
{
	stk::io::PointSetOutputStream<DIM, POS, VAL> ptsStream;
	ptsStream.setWeight(i_weight);
	ptsStream.open(i_filename, true, true);
	for(unsigned int i=0; i<o_ptsList.size(); i++)
	{
		ptsStream.write(o_ptsList.at(i));
	}
}

template<int DIM, typename POS, typename VAL>
void write(const std::string i_filename, const std::vector< PointSet<DIM, POS, VAL> >& o_ptsList, const bool& i_weight = true)
{ writePts<DIM, POS, VAL>(i_filename, o_ptsList, i_weight); }

template<int DIM, typename POS, typename VAL>
void writePts(const std::string i_filename, const PointSet<DIM, POS, VAL>& o_pts, const bool& i_weight = true)
{
	stk::io::PointSetOutputStream<DIM, POS, VAL> ptsStream;
	ptsStream.setWeight(i_weight);
	ptsStream.open(i_filename, false, true);
	ptsStream.write(o_pts);
}

template<int DIM, typename POS, typename VAL>
void write(const std::string i_filename, const PointSet<DIM, POS, VAL>& o_pts, const bool& i_weight = true)
{ writePts<DIM, POS, VAL>(i_filename, o_pts, i_weight); }

template<int DIM, typename POS, typename VAL>
void writeDat(const std::string i_filename, const std::vector< PointSet<DIM, POS, VAL> >& o_ptsList, const bool& i_weight = false)
{
	stk::io::PointSetOutputStream<DIM, POS, VAL> ptsStream;
	ptsStream.setWeight(i_weight);
	ptsStream.open(i_filename, true, false);
	for(unsigned int i=0; i<o_ptsList.size(); i++)
	{
		ptsStream.write(o_ptsList.at(i));
	}
}

template<int DIM, typename POS, typename VAL>
void writeDat(const std::string i_filename, const PointSet<DIM, POS, VAL>& o_pts, const bool& i_weight = false)
{
	stk::io::PointSetOutputStream<DIM, POS, VAL> ptsStream;
	ptsStream.setWeight(i_weight);
	ptsStream.open(i_filename, false, false);
	ptsStream.write(o_pts);
}



#ifdef CAIRO_ENABLED

template<typename VAL = double>
class PointSetGraphics
{
	public:
		enum PointSetGraphicsFlag
		{
			Default = 0,
			Wrong,
			Fixed,
			Hot,
			Cold
		};
		enum PointSetGraphicsLabel
		{
			None,
			Id,
			Value
		};
		enum PointSetGraphicsVoronoiColor
		{
			VOR_NONE,
			VOR_CONNEXITY,
			VOR_AREA,
			VOR_CLASS_BOUND,
			VOR_CLASS_COLOR
		};
	
	protected:
		std::map<float, const PointSet<2, double, VAL>*> m_timeTrace;
		std::map<unsigned int, PointSetGraphicsFlag> m_flags;
		const PointSet<2, double, VAL>* m_pointset;
		int m_borderSize;
		bool m_computeDelaunay;
		bool m_delaunay;
		bool m_voronoi;
		bool m_timeTracePoint;
		PointSetGraphicsLabel m_labelType;
		PointSetGraphicsVoronoiColor m_voronoiColorType;
		bool m_zoom;
		Vector2d m_zoomCenter;
		double m_zoomScale;
		double m_pointRadius;
		
		class DomainToImage
		{
			protected:
				Vector2d imgSize;
				
			public:
				DomainToImage(const Vector2d& size);
				virtual Vector2d convert(const Vector2d& p) const = 0;
		};
		
		class DomainToImageGlobal : public DomainToImage
		{
			protected:
				Vector2d viewMin;
				Vector2d viewSize;
				
			public:
				DomainToImageGlobal(const Vector2d& size, const Domain<2, double>* domain, int borderSize);
				virtual Vector2d convert(const Vector2d& p) const;
		};
		
		class DomainToImageZoom : public DomainToImage
		{
			protected:
				double scale;
				Vector2d center;
				Vector2d viewMin;
				Vector2d viewSize;
				
			public:
				DomainToImageZoom(const Vector2d& size, const Domain<2, double>* domain, int borderSize, const Vector2d& c, double s);
				virtual Vector2d convert(const Vector2d& p) const;
		};
		
	public:
		PointSetGraphics()
		{
			m_borderSize = 26;
			m_voronoi = false;
			m_delaunay = false;
			m_timeTracePoint = false;
			m_labelType = None;
			m_zoom = false;
			m_zoomScale = 1.0;
			m_computeDelaunay = false;
			m_pointRadius = 1.0;
		}
		
		void radiusSize(double i_radius)
		{
			m_pointRadius = i_radius;
		}
		
		void borderSize(int i_size)
		{
			m_borderSize = i_size;
		}
		
		int borderSize() const
		{
			return m_borderSize;
		}
		
		void timeTrace(const PointSet<2, double, VAL>& i_pointset, float i_time)
		{
			m_timeTrace[i_time] = &i_pointset;
		}
		
		void pointset(const PointSet<2, double, VAL>& i_pointset)
		{
			m_pointset = &i_pointset;
		}
		
		void addFlag(unsigned int id, const PointSetGraphicsFlag& flag)
		{
			m_flags[id] = flag;
		}
		
		void enableDelaunay()
		{
			m_delaunay = true;
			m_computeDelaunay = true;
		}
		
		void enableVoronoi(const PointSetGraphicsVoronoiColor& voronoiColor = VOR_NONE)
		{
			m_voronoi = true;
			m_voronoiColorType = voronoiColor;
			m_computeDelaunay = true;
		}
		
		void enableTimeTracePoint()
		{
			m_timeTracePoint = true;
		}
		
		void enableLabel(const PointSetGraphicsLabel& labelType)
		{
			m_labelType = labelType;
		}
		
		void enableZoom(const Vector2d& center, const double& size)
		{
			m_zoomCenter = center;
			m_zoom = true;
			m_zoomScale = size;
		}
		
		void draw(const std::string& i_filename, const Vector2i& i_size) const;
};

void draw(std::string i_filename, const PointSet2dd& i_pts, const Vector2i& i_size);
void draw(std::string i_filename, const PointSet2di& i_pts, const Vector2i& i_size);

void draw(
	std::string i_filename,
	const PointSet2dd& i_pts,
	const Vector2i& i_size);

void draw(
	std::string i_filename,
	const PointSet2dc& i_pts,
	const Vector2i& i_size);

void draw(
	std::string i_filename,
	const PointSet2di& i_pts,
	const Vector2i& i_size);

void drawMotion(
	std::string i_filename,
	const PointSet2dd& i_before,
	const PointSet2dd& i_after,
	const Vector2i& i_size);
	
void drawWeight(
	std::string i_filename,
	const PointSet2dd& i_pts,
	const Vector2i& i_size);
#endif

}

}

#endif
