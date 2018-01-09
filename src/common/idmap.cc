/*
 * FUNC - Functional Analysis of Gene Expression Data
 * Copyright (C) 2002  Bjoern Muetzel, Kay Pruefer
 * 
 * This program is modifiable/redistributable under the terms of the
 * GNU General Public License.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; see the file COPYING. If not, write to the
 * Free Software Foundation, Inc., 59 Temple Place - Suite 330,
 * Boston, MA 02111-1307, USA.
 */




#include "idmap.h"
#include <sstream>

idmap::idmap( istream &in ) 
{
	string line ;

	while ( in ) {
		// first row: id
		getline( in, line, '\t' ) ;
		string id( line ) ;
		
		// skip 2 fields
		getline( in, line, '\t' ) ;
		getline( in, line, '\t' ) ;

		// GO:Number
		getline( in, line, '\t' ) ;
		string go( line ) ;
		
		if ( id.size() > 0 && go.size() > 0 ) (*this)[id]=go ;
		
		getline( in, line, '\n' ) ;
	}
}

string idmap::get_id_for_go( string &go ) 
{
	for ( map<string,string>::const_iterator it = this->begin() ; 
			it != this->end() ; ++it ) 
			if ( it->second == go ) return it->first ;

	return 0 ;
}
