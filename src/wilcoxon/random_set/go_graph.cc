
#include "go_graph.h"
#include <sstream>

using std::cerr;
using std::endl;
using std::istringstream;

go_graph::go_graph( set<string> &nodes, istream &term2term, idmap &idm_ )
	: idm( idm_ ) 
{
	// read nodes, make go_objs
	map<string, go_obj*> temp_graph ;
	for ( set<string>::const_iterator it = nodes.begin() ; it!=nodes.end() ;
			++it ) 
	{
		temp_graph[*it] = new go_obj( idm[*it] ) ;
	}

	// read term2term file, create graph
	while ( term2term ) {
		// 4 columns, 1 a id, 2 type of relationship, 3 parent, 4 child
		// new format: 5 col. Last column: 5 "complete"
		string s1 ;

		// skip first 2 columns
		getline( term2term, s1, '\t' ) ;
		getline( term2term, s1, '\t' ) ;


		// parent id
		getline( term2term, s1, '\t' ) ;

		/* The format of the term_db_tables has changed in 12.2004. 
		   a new column appeared in the term2term table. So we have to check
		   whether we are using the old or the new format und extract
		   the child_id */

		map<string, go_obj*>::const_iterator par = temp_graph.find( s1 ) ;
		if ( par != temp_graph.end() ) {
			// child id
			getline( term2term, s1, '\n' ) ;
			string child_id ;
			string::size_type tab_pos = s1.find( '\t' ) ;
			if ( tab_pos != string::npos ) {
				child_id = s1.substr( 0, tab_pos ) ;
			} else child_id = s1 ;
			map<string, go_obj*>::const_iterator child = 
						temp_graph.find( child_id ) ;
			if ( child != temp_graph.end() ) {
				child->second->add_parent( par->second ) ;
			//	par->second->add_child( child->second ) ;
			}	
		} else {
			// skip rest if the parent node is not part of the graph
			getline( term2term, s1, '\n' ) ;
		}
	}
	// rewrite map file because detectedfile has
	// go_ids instead of database ids.
	for ( map<string, go_obj*>::const_iterator it = temp_graph.begin() ;
			it != temp_graph.end() ; ++it ) 
	{
		graph[idm[it->first]] = it->second ;
	}
}

go_graph::~go_graph(  )
{
	for ( map<string, go_obj*>::const_iterator it = graph.begin() ;
			it != graph.end() ; ++it ) 
	{
		delete it->second ;
	}
}

void go_graph::get_parents( string &go, set<go_obj*>* parents ) {
	
	if ( graph.find( go ) != graph.end() ) {
		graph[go]->get_parents( parents ) ;
	} else {
		cerr << "Cannot find " << go << endl ;
	} 
}


void go_graph::clear_genes(  ) {
	for ( map<string, go_obj*>::const_iterator it = graph.begin() ;
			it != graph.end() ; ++it ) 
	{
		it->second->clear_genes( ) ;
	}
}

void go_graph::print_header( ostream &os ) 
{
	for ( map<string,go_obj*>::const_iterator it = graph.begin() ;
				it != graph.end() ; ++it ) 
	{
		os << it->first << '\t' ;
	}
	os << '\n' ;
	for ( map<string,go_obj*>::const_iterator it = graph.begin() ;
				it != graph.end() ; ++it ) 
	{
		it->second->print_n( os ) ;
	}
	os << '\n' ;
}

void go_graph::print_sumranks( ostream &os )
{
	for ( map<string,go_obj*>::const_iterator it = graph.begin() ;
				it != graph.end() ; ++it ) 
	{
		it->second->print_sumranks( os ) ; 
	}
	os << '\n' ;
}
