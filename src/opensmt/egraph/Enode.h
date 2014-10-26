/*********************************************************************
Author: Roberto Bruttomesso <roberto.bruttomesso@gmail.com>

OpenSMT -- Copyright (C) 2010, Roberto Bruttomesso

OpenSMT is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

OpenSMT is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with OpenSMT. If not, see <http://www.gnu.org/licenses/>.
*********************************************************************/

#ifndef ENODE_H
#define ENODE_H

#include "egraph/EnodeTypes.h"
#include "common/Otl.h"
#include <limits>
#include <unordered_set>

class Enode
{
public:

  //
  // Constructor for Enil
  //
  Enode  ( );
  //
  // Constructor for symbols
  //
  Enode ( const enodeid_t       // id
        , const char *          // name/value
        , const etype_t         // enode type
        , Snode *               // Sort
        );
  //
  // Constructor for terms and lists
  //
  Enode ( const enodeid_t       // id
        , Enode *               // car
        , Enode *               // cdr
        );
  //
  // Constructor for defs
  //
  Enode ( const enodeid_t       // id
        , Enode *               // def
        );
  //
  // Destructor
  //
  ~Enode ( );
  //
  // Check if a node is Enil
  //
  inline bool isEnil            ( ) const { return id == ENODE_ID_ENIL; }
  inline bool isList            ( ) const { return (properties & ETYPE_MASK) == ETYPE_LIST; }
  inline bool isTerm            ( ) const { return (properties & ETYPE_MASK) == ETYPE_TERM; }
  inline bool isSymb            ( ) const { return (properties & ETYPE_MASK) == ETYPE_SYMB; }
  inline bool isNumb            ( ) const { return (properties & ETYPE_MASK) == ETYPE_NUMB; }
  inline bool isDef             ( ) const { return (properties & ETYPE_MASK) == ETYPE_DEF; }

  inline void setEtype          ( const etype_t t )
  {
    assert( t == ETYPE_SYMB
         || t == ETYPE_NUMB
         || t == ETYPE_LIST
         || t == ETYPE_TERM
         || t == ETYPE_DEF );
    properties |= t;
  }

  inline void setArity            ( const unsigned a ) { assert( a <= ARITY_N ); properties |= (a << ARITY_SHIFT); }
  //
  // Check if a term node represents a certain symbol
  //
  inline bool isPlus              ( ) const { return hasSymbolId( ENODE_ID_PLUS        ); }
  inline bool isMinus             ( ) const { return hasSymbolId( ENODE_ID_MINUS       ); }
  inline bool isUminus            ( ) const { return hasSymbolId( ENODE_ID_UMINUS      ); }
  inline bool isTimes             ( ) const { return hasSymbolId( ENODE_ID_TIMES       ); }
  inline bool isDiv               ( ) const { return hasSymbolId( ENODE_ID_DIV         ); }

  /* added for dReal2 */
  inline bool isAcos              ( ) const { return hasSymbolId( ENODE_ID_ACOS              ); }
  inline bool isAsin              ( ) const { return hasSymbolId( ENODE_ID_ASIN              ); }
  inline bool isAtan              ( ) const { return hasSymbolId( ENODE_ID_ATAN              ); }
  inline bool isAtan2             ( ) const { return hasSymbolId( ENODE_ID_ATAN2             ); }
  inline bool isMatan             ( ) const { return hasSymbolId( ENODE_ID_MATAN             ); }
  inline bool isSqrt 	          ( ) const { return hasSymbolId( ENODE_ID_SQRT) ;}
  inline bool isSafeSqrt          ( ) const { return hasSymbolId( ENODE_ID_SAFESQRT      ); }
  inline bool isExp               ( ) const { return hasSymbolId( ENODE_ID_EXP         ); }
  inline bool isLog               ( ) const { return hasSymbolId( ENODE_ID_LOG         ); }
  inline bool isPow               ( ) const { return hasSymbolId( ENODE_ID_POW         ); }
  inline bool isAbs               ( ) const { return hasSymbolId( ENODE_ID_ABS         ); }
  inline bool isSin               ( ) const { return hasSymbolId( ENODE_ID_SIN         ); }
  inline bool isCos               ( ) const { return hasSymbolId( ENODE_ID_COS         ); }
  inline bool isTan               ( ) const { return hasSymbolId( ENODE_ID_TAN         ); }
  inline bool isSinh              ( ) const { return hasSymbolId( ENODE_ID_SINH        ); }
  inline bool isCosh              ( ) const { return hasSymbolId( ENODE_ID_COSH        ); }
  inline bool isTanh              ( ) const { return hasSymbolId( ENODE_ID_TANH        ); }
  /* ---------------------- */

  inline bool isEq                ( ) const { return hasSymbolId( ENODE_ID_EQ          ); }
  inline bool isLeq               ( ) const { return hasSymbolId( ENODE_ID_LEQ         ); }
  inline bool isGeq               ( ) const { return hasSymbolId( ENODE_ID_GEQ         ); }
  inline bool isLt                ( ) const { return hasSymbolId( ENODE_ID_LT          ); }
  inline bool isGt                ( ) const { return hasSymbolId( ENODE_ID_GT          ); }
  inline bool isStore             ( ) const { return hasSymbolId( ENODE_ID_STORE       ); }
  inline bool isSelect            ( ) const { return hasSymbolId( ENODE_ID_SELECT      ); }
  inline bool isImplies           ( ) const { return hasSymbolId( ENODE_ID_IMPLIES     ); }
  inline bool isAnd               ( ) const { return hasSymbolId( ENODE_ID_AND         ); }
  inline bool isOr                ( ) const { return hasSymbolId( ENODE_ID_OR          ); }
  inline bool isNot               ( ) const { return hasSymbolId( ENODE_ID_NOT         ); }
  inline bool isXor               ( ) const { return hasSymbolId( ENODE_ID_XOR         ); }
  inline bool isIff               ( ) const { return hasSymbolId( ENODE_ID_EQ ) && get1st( )->hasSortBool( ); }
  inline bool isTrue              ( ) const { return hasSymbolId( ENODE_ID_TRUE        ); }
  inline bool isFalse             ( ) const { return hasSymbolId( ENODE_ID_FALSE       ); }
  inline bool isIte               ( ) const { return hasSymbolId( ENODE_ID_ITE         ); }
  inline bool isDistinct          ( ) const { return hasSymbolId( ENODE_ID_DISTINCT    ); }
  inline bool isForallT           ( ) const { return hasSymbolId( ENODE_ID_FORALLT     ); }
  inline bool isIntegral          ( ) const { return hasSymbolId( ENODE_ID_INTEGRAL    ); }
  inline bool isConnect		        ( ) const { return hasSymbolId( ENODE_ID_CONNECT     ); }
  inline bool isPIntegral         ( ) const { return hasSymbolId( ENODE_ID_PINTEGRAL    ); }
  /*
  inline bool isBvslt             ( ) const { return hasSymbolId( ENODE_ID_BVSLT       ); }
  inline bool isBvsgt             ( ) const { return hasSymbolId( ENODE_ID_BVSGT       ); }
  inline bool isBvsle             ( ) const { return hasSymbolId( ENODE_ID_BVSLE       ); }
  inline bool isBvsge             ( ) const { return hasSymbolId( ENODE_ID_BVSGE       ); }
  inline bool isBvult             ( ) const { return hasSymbolId( ENODE_ID_BVULT       ); }
  inline bool isBvugt             ( ) const { return hasSymbolId( ENODE_ID_BVUGT       ); }
  inline bool isBvule             ( ) const { return hasSymbolId( ENODE_ID_BVULE       ); }
  inline bool isBvuge             ( ) const { return hasSymbolId( ENODE_ID_BVUGE       ); }
  inline bool isConcat            ( ) const { return hasSymbolId( ENODE_ID_CONCAT      ); }
  inline bool isCbe               ( ) const { return hasSymbolId( ENODE_ID_CBE         ); }
  inline bool isBvand             ( ) const { return hasSymbolId( ENODE_ID_BVAND       ); }
  inline bool isBvor              ( ) const { return hasSymbolId( ENODE_ID_BVOR        ); }
  inline bool isBvxor             ( ) const { return hasSymbolId( ENODE_ID_BVXOR       ); }
  inline bool isBvnot             ( ) const { return hasSymbolId( ENODE_ID_BVNOT       ); }
  inline bool isBvadd             ( ) const { return hasSymbolId( ENODE_ID_BVADD       ); }
  inline bool isBvsub             ( ) const { return hasSymbolId( ENODE_ID_BVSUB       ); }
  inline bool isBvmul             ( ) const { return hasSymbolId( ENODE_ID_BVMUL       ); }
  inline bool isBvneg             ( ) const { return hasSymbolId( ENODE_ID_BVNEG       ); }
  inline bool isBvlshr            ( ) const { return hasSymbolId( ENODE_ID_BVLSHR      ); }
  inline bool isBvashr            ( ) const { return hasSymbolId( ENODE_ID_BVASHR      ); }
  inline bool isBvshl             ( ) const { return hasSymbolId( ENODE_ID_BVSHL       ); }
  inline bool isBvsrem            ( ) const { return hasSymbolId( ENODE_ID_BVSREM      ); }
  inline bool isBvurem            ( ) const { return hasSymbolId( ENODE_ID_BVUREM      ); }
  inline bool isBvsdiv            ( ) const { return hasSymbolId( ENODE_ID_BVSDIV      ); }
  inline bool isBvudiv            ( ) const { return hasSymbolId( ENODE_ID_BVUDIV      ); }
  inline bool isSignExtend        ( ) const { int i; return sscanf( car->getName( ), "sign_extend[%d]", &i ) == 1; }
  bool        isSignExtend        ( int * );
  inline bool isZeroExtend        ( ) const { return hasSymbolId( ENODE_ID_ZERO_EXTEND ); }
  inline bool isBoolcast          ( ) const { return hasSymbolId( ENODE_ID_BOOLCAST    ); }
  inline bool isWord1cast         ( ) const { return hasSymbolId( ENODE_ID_WORD1CAST   ); }
  */
  inline bool isUp                ( ) const { return car->id > ENODE_ID_LAST && isAtom( ) && getArity( ) > 0; }
  inline bool isUf                ( ) const { return car->id > ENODE_ID_LAST && !isAtom( ) && getArity( ) > 0; }

  //  inline bool isCostIncur       ( ) const { return hasSymbolId( ENODE_ID_CTINCUR ); }
  //  inline bool isCostBound       ( ) const { return hasSymbolId( ENODE_ID_CTBOUND ); }

  bool        isVar               ( ) const; // True if it is a variable
  bool        isConstant          ( ) const; // True if it is a constant
  bool        isLit               ( ) const; // True if it is a literal
  bool        isAtom              ( ) const; // True if it is an atom
  bool        isTLit              ( ) const; // True if it is a theory literal
  bool        isTAtom             ( ) const; // True if it is a theory atom
  bool        isBooleanOperator   ( ) const; // True if it is a boolean operator
  bool        isArithmeticOp      ( ) const; // True if root is an arith term
  bool        isUFOp              ( ) const; // True if root is UF

  inline bool hasSortBool( ) const
  {
    assert( isTerm( ) );
    if ( isIte( ) ) return get2nd( )->hasSortBool( );
    return car->symb_data->sort->hasSortBool( );
  }
  inline bool hasSortReal( ) const
  {
    assert( isTerm( ) );
    if ( isIte( ) ) return get2nd( )->hasSortReal( );
    return car->symb_data->sort->hasSortReal( );
  }
  inline bool hasSortInt( ) const
  {
    assert( isTerm( ) );
    if ( isIte( ) ) return get2nd( )->hasSortInt( );
    return car->symb_data->sort->hasSortInt( );
  }
  inline bool hasSortArray( ) const
  {
    assert( isTerm( ) );
    if ( isIte( ) ) return get2nd( )->hasSortArray( );
    return car->symb_data->sort->hasSortArray( );
  }
  inline bool hasSortIndex( ) const
  {
    assert( isTerm( ) );
    if ( isIte( ) ) return get2nd( )->hasSortIndex( );
    return car->symb_data->sort->hasSortIndex( );
  }
  inline bool hasSortElem( ) const
  {
    assert( isTerm( ) );
    if ( isIte( ) ) return get2nd( )->hasSortElem( );
    return car->symb_data->sort->hasSortElem( );
  }
  inline bool hasSortUndef( ) const
  {
    assert( isTerm( ) );
    if ( isIte( ) ) return get2nd( )->hasSortUndef( );
    return car->symb_data->sort->hasSortUndef( );
  }

  inline bool hasCongData         ( ) const { return cong_data != NULL; }

  void        allocCongData       ( );
  void        deallocCongData     ( );

  bool 		isHolder;

  //
  // Getty and Setty methods
  //
  inline enodeid_t            getId      ( ) const { return id; }
  inline unsigned             getArity   ( ) const { return ((properties & ARITY_MASK) >> ARITY_SHIFT); }
  Snode *                     getSort    ( ) const { assert( isTerm( ) || isSymb( ) ); return isTerm( ) ? car->symb_data->sort : symb_data->sort; }
  Snode *                     getLastSort( );
  inline string   getName                ( ) const { assert( isSymb( ) || isNumb( ) ); assert( symb_data ); return stripName( symb_data->name ); }
  inline string   getNameFull            ( ) const { assert( isSymb( ) || isNumb( ) ); assert( symb_data ); return symb_data->name; }
  inline Enode *  getCar                 ( ) const { return car; }
  inline Enode *  getCdr                 ( ) const { return cdr; }
  inline Enode *  getDef                 ( ) const { assert( isDef( ) ); assert( car ); return car; }

  inline Enode *  getNext                ( ) const { assert( isTerm( ) || isList( ) ); assert( cong_data ); return cong_data->next; }
  inline int      getSize                ( ) const { assert( isTerm( ) || isList( ) ); assert( cong_data ); return cong_data->size; }
  inline Enode *  getParent              ( ) const { assert( isTerm( ) || isList( ) ); assert( cong_data ); return cong_data->parent; }
  inline Enode *  getSameCar             ( ) const { assert( isTerm( ) || isList( ) ); assert( cong_data ); return cong_data->same_car; }
  inline Enode *  getSameCdr             ( ) const { assert( isTerm( ) || isList( ) ); assert( cong_data ); return cong_data->same_cdr; }
  inline int      getParentSize          ( ) const { assert( isTerm( ) || isList( ) ); assert( cong_data ); return cong_data->parent_size; }
  inline Enode *  getCgPtr               ( ) const { assert( isTerm( ) || isList( ) ); assert( cong_data ); return cong_data->cg_ptr; }
  inline Elist *  getForbid              ( ) const { assert( isTerm( ) || isList( ) ); assert( cong_data ); return cong_data->forbid; }
  inline dist_t   getDistClasses         ( ) const { assert( isTerm( ) || isList( ) ); assert( cong_data ); return cong_data->dist_classes; }

  double          getValue               ( ) const;
  std::unordered_set<Enode *> get_vars   ( );
  std::unordered_set<Enode *> get_constants   ( );

  double          getLowerBound          ( ) const; //added for dReal2
  double          getUpperBound          ( ) const; //added for dReal2
  double          getPrecision           ( ) const; //added for dReal2
  bool            hasPrecision           ( ) const; //added for dReal2
  double          getComplexValue        ( ) const;
  void            setValue               ( const double );
  void            setLowerBound          ( const double ); //added for dReal2
  void            setUpperBound          ( const double ); //added for dReal2
  void            setPrecision           ( const double ); //added for dReal2
  bool            hasValue               ( ) const;
  Enode *         getRoot                ( ) const;
  enodeid_t       getCid                 ( ) const;
  Enode *         getConstant            ( ) const;

  inline Enode *  getExpReason           ( ) const { assert( isTerm( ) && cong_data && cong_data->term_data ); return cong_data->term_data->exp_reason; }
  inline Enode *  getExpParent           ( ) const { assert( isTerm( ) && cong_data && cong_data->term_data ); return cong_data->term_data->exp_parent; }
  inline Enode *  getExpRoot             ( ) const { assert( isTerm( ) && cong_data && cong_data->term_data ); return cong_data->term_data->exp_root; }
  inline int      getExpClassSize        ( ) const { assert( isTerm( ) && cong_data && cong_data->term_data ); return cong_data->term_data->exp_class_size; }
  inline Enode *  getExpHighestNode      ( ) const { assert( isTerm( ) && cong_data && cong_data->term_data ); return cong_data->term_data->exp_highest_node; }
  // inline Reason * getExpReason           ( ) const { assert( isTerm( ) && cong_data && cong_data->term_data ); return cong_data->term_data->exp_reason; }
  inline int      getExpTimeStamp        ( ) const { assert( isTerm( ) && cong_data && cong_data->term_data ); return cong_data->term_data->exp_time_stamp; }

  inline lbool    getPolarity            ( ) const { assert( isTerm( ) && atom_data ); return atom_data->polarity; }
  inline bool     hasPolarity            ( ) const { assert( isTerm( ) && atom_data ); return atom_data->has_polarity; }
  inline lbool    getDeduced             ( ) const { assert( isTerm( ) && atom_data ); return atom_data->deduced; }
  inline bool     isDeduced              ( ) const { assert( isTerm( ) && atom_data ); return atom_data->is_deduced; }
  inline lbool    getDecPolarity         ( )       { assert( isAtom( ) && atom_data ); return atom_data->dec_polarity; }
  inline int      getWeightInc           ( )       { assert( isAtom( ) && atom_data ); return atom_data->weight_inc; }
  inline int      getDedIndex            ( ) const { assert( isTerm( ) && atom_data ); return atom_data->ded_index; }
  inline int      getDistIndex           ( ) const { assert( isTerm( ) && atom_data ); return atom_data->dist_index; }

  inline Enode *  getCb                  ( ) const { assert( isTerm( ) && cong_data && cong_data->term_data ); return cong_data->term_data->cb; }
  inline Enode *  getRef                 ( ) const { assert( isTerm( ) && cong_data && cong_data->term_data ); return cong_data->root; }
  inline int      getWidth               ( ) const { assert( isTerm( ) || isSymb( ) || isNumb( ) ); return (properties & WIDTH_MASK); }

  inline void     setWidth               ( const uint32_t w )
  {
    assert( isTerm( ) );
    assert( w < MAX_WIDTH );
    // Reset width
    properties &= ~WIDTH_MASK;
    properties |= w;
    assert( getWidth( ) == static_cast<int>( w ) );
  }

  inline void    setRoot                ( Enode * e )        { assert( isTerm( ) || isList( ) ); assert( cong_data ); cong_data->root = e; }
  inline void    setCid                 ( const enodeid_t c ){ assert( isTerm( ) || isList( ) ); assert( cong_data ); cong_data->cid = c; }
  inline void    setDef                 ( Enode * e )        { assert( e ); assert( isDef( ) ); car = e; }
  inline void    setNext                ( Enode * e )        { assert( isTerm( ) || isList( ) ); assert( cong_data ); cong_data->next = e; }
  inline void    setSize                ( const int s )      { assert( isTerm( ) || isList( ) ); assert( cong_data ); cong_data->size = s; }
  inline void    setParent              ( Enode * e )        { assert( isTerm( ) || isList( ) ); assert( cong_data ); cong_data->parent = e; }
  inline void    setSameCar             ( Enode * e )        { assert( isTerm( ) || isList( ) ); assert( cong_data ); cong_data->same_car = e; }
  inline void    setSameCdr             ( Enode * e )        { assert( isTerm( ) || isList( ) ); assert( cong_data ); cong_data->same_cdr = e; }
  inline void    setParentSize          ( const int s )      { assert( isTerm( ) || isList( ) ); assert( cong_data ); cong_data->parent_size = s; }
  inline void    setCgPtr               ( Enode * e )        { assert( isTerm( ) || isList( ) ); assert( cong_data ); cong_data->cg_ptr = e; }
  inline void    setForbid              ( Elist * l )        { assert( isTerm( ) || isList( ) ); assert( cong_data ); cong_data->forbid = l; }
  inline void    setDistClasses         ( const dist_t & d ) { assert( isTerm( ) || isList( ) ); assert( cong_data ); cong_data->dist_classes = d; }

  inline void    setConstant            ( Enode * e )            { assert( isTerm( ) && cong_data && cong_data->term_data ); assert( e == NULL || e->isConstant( ) ); cong_data->term_data->constant = e; }
  inline void    setExpReason           ( Enode * r )            { assert( isTerm( ) && cong_data && cong_data->term_data ); cong_data->term_data->exp_reason = r; }
  inline void    setExpParent           ( Enode * e )            { assert( isTerm( ) && cong_data && cong_data->term_data ); cong_data->term_data->exp_parent = e; }
  inline void    setExpRoot             ( Enode * e )            { assert( isTerm( ) && cong_data && cong_data->term_data ); cong_data->term_data->exp_root = e; }
  inline void    setExpClassSize        ( const int s )          { assert( isTerm( ) && cong_data && cong_data->term_data ); cong_data->term_data->exp_class_size = s; }
  inline void    setExpHighestNode      ( Enode * e )            { assert( isTerm( ) && cong_data && cong_data->term_data ); cong_data->term_data->exp_highest_node = e; }
  inline void    setExpTimeStamp        ( const int t )          { assert( isTerm( ) && cong_data && cong_data->term_data ); cong_data->term_data->exp_time_stamp = t; }
  inline void    setPolarity            ( const lbool p )        { assert( isTerm( ) && atom_data ); assert( !atom_data->has_polarity ); atom_data->polarity = p; atom_data->has_polarity = true; }
  inline void    resetPolarity          ( )                      { assert( isTerm( ) && atom_data ); assert( atom_data->has_polarity ); atom_data->has_polarity = false; }
  inline void    setDeduced             ( const lbool d, int i ) { assert( isTerm( ) && atom_data ); assert( !atom_data->is_deduced ); atom_data->deduced = d; atom_data->ded_index = i; atom_data->is_deduced = true; }
  inline void    setDeduced             ( const bool s, int i )  { setDeduced( s ? l_False : l_True, i ); }
  inline void    resetDeduced           ( )                      { assert( isTerm( ) && atom_data ); assert( atom_data->is_deduced ); atom_data->is_deduced = false; }
  inline void    setDecPolarity         ( const lbool s )        { assert( isAtom( ) && atom_data ); atom_data->dec_polarity = s; }
  inline void    setWeightInc           ( const int w )          { assert( isAtom( ) && atom_data ); atom_data->weight_inc = w; }
  inline void    setDistIndex           ( const int d )          { assert( isTerm( ) && atom_data ); atom_data->dist_index = d; }

  inline void    setCb                  ( Enode * e )            { assert( isTerm( ) && cong_data && cong_data->term_data ); cong_data->term_data->cb = e; }

#if 0
  inline void    setDynamic             ( Enode * e )            { assert( dynamic == NULL ); dynamic = e; }
  inline void    resetDynamic           ( )                      { assert( dynamic != NULL ); dynamic = NULL; }
  inline Enode * getDynamic             ( ) const                { return dynamic; }
#endif

  //
  // Congruence closure related routines
  //
  void           addParent              ( Enode * );   // Adds a parent to this node. Useful for congruence closure
  void           addParentTail          ( Enode * );   // Adds a parent to the tail of this node. Useful for congruence closure
  void           removeParent           ( Enode * );   // Remove a parent of this node
  //
  // Shortcuts for retrieving a term's arguments
  //
  inline Enode * get1st                 ( ) const;     // Get first argument in constant time
  inline Enode * get2nd                 ( ) const;     // Get second argument in constant time
  inline Enode * get3rd                 ( ) const;     // Get third argument in constant time

  bool           addToCongruence        ( ) const;
  unsigned       sizeInMem              ( ) const;

  void           print_infix( ostream & os, lbool polarity, string const & variable_postfix = "") const;
  void           print                  ( ostream & ); // Prints the

  string         stripName              ( string ) const;
  void           printSig               ( ostream & ); // Prints the enode signature

#ifdef BUILD_64
  inline enodeid_pair_t          getSig    ( ) const { return encode( car->getRoot( )->getCid( ), cdr->getRoot( )->getCid( ) ); }
#else
  inline const Pair( enodeid_t ) getSig    ( ) const { return make_pair( car->getRoot( )->getCid( ), cdr->getRoot( )->getCid( ) ); }
#endif
  inline enodeid_t               getSigCar ( ) const { return car->getRoot( )->getCid( ); }
  inline enodeid_t               getSigCdr ( ) const { return cdr->getRoot( )->getCid( ); }

  inline friend ostream &  operator<<( ostream & os, Enode * e )    { assert( e ); e->print( os ); return os; }

  struct idLessThan
  {
    inline bool operator( )( Enode * x, Enode * y )
    {
      return (x->getCar( )->getId( ) <  y->getCar( )->getId( ))
          || (x->getCar( )->getId( ) == y->getCar( )->getId( ) && x->getCdr( )->getId( ) < y->getCdr( )->getId( ) );
    }
  };

  struct cidLessThan
  {
    inline bool operator( )( Enode * x, Enode * y )
    {
      if ( x == y ) return false;
      if ( x->isEnil( ) ) return true;
      if ( y->isEnil( ) ) return false;
      return (x->getCar( )->getCid( ) <  y->getCar( )->getCid( ))
          || (x->getCar( )->getCid( ) == y->getCar( )->getCid( ) && x->getCdr( )->getCid( ) < y->getCdr( )->getCid( ) );
    }
  };

private:
  //
  // Standard informations for terms
  //
  const enodeid_t   id;         // Node unique identifier
  uint32_t          properties; // Contains all the properties of this node (see EnodeTypes.h for bitfields definition)
  Enode *           car;        // For car / defs
  Enode *           cdr;        // For cdr
  union {
  CongData *        cong_data;  // For terms and lists that go in the congruence
  SymbData *        symb_data;  // For symbols/numbers
  };
  AtomData *        atom_data;  // For atom terms only
  double *          value;      // enode value (modified for dReal2)
  double *          lb;         // enode lower bound
  double *          ub;         // enode upper bound
  double            precision; //added for dReal2

#if 0
  Enode *           dynamic;    // Pointer to dynamic equivalent
#endif

  inline bool       hasSymbolId    ( const enodeid_t id ) const { assert( isTerm( ) ); return car->getId( ) == id; }
};

inline double       Enode::getValue ( ) const
{
  assert( hasValue( ) );
  return *value;
}

inline double Enode::getLowerBound ( ) const
{
    if (lb != NULL)
        return *lb;
    else
        return -std::numeric_limits<double>::infinity();
}

inline double Enode::getUpperBound ( ) const
{
    if (ub != NULL)
        return *ub;
    else
        return +std::numeric_limits<double>::infinity();
}

inline double Enode::getPrecision ( ) const
{
  assert( hasPrecision() );
  return precision;
}

inline bool Enode::hasPrecision ( ) const
{
  return precision != 0.0;
}



inline double Enode::getComplexValue( ) const
{
  if( isDiv( ) )
    return get1st( )->getCar( )->getComplexValue( ) / get2nd( )->getCar( )->getComplexValue( );
  else if( isUminus( ) )
    return -(*value);
  else
    return *value;
}

inline void Enode::setValue ( const double v )
{
  assert( isTerm( ) );
  *value = v;
}

inline void Enode::setLowerBound ( const double v )
{
  assert( isTerm( ) );
  lb = new double;
  *lb = v;
}

inline void Enode::setUpperBound ( const double v )
{
  assert( isTerm( ) );
  ub = new double;
  *ub = v;
}

inline void Enode::setPrecision ( const double v )
{
  assert( isTerm( ) );

  if( isNot() )
  {
    //only set precision on the atom (not the literal)
    get1st()->setPrecision(v);
  }
  else {
    precision = v;
  }
}

inline bool Enode::hasValue( ) const
{
  assert( isTerm( ) );
  return value != NULL;
}

inline Enode * Enode::getRoot ( ) const
{
  assert( !isDef( ) );
  if ( (isTerm( ) || isList( )) && cong_data != NULL )
    return cong_data->root;
  return const_cast<Enode *>(this);
}

inline enodeid_t Enode::getCid ( ) const
{
  assert( !isDef( ) );
  if ( (isTerm( ) || isList( )) && cong_data )
    return cong_data->cid;
  return id;
}

inline Enode * Enode::getConstant ( ) const
{
  assert( isTerm( ) || isList( ) );
  // assert( cong_data );
  if ( cong_data == NULL )
    return NULL;

  if ( isTerm( ) )
  {
    assert( cong_data );
    assert( cong_data->term_data );
    return cong_data->term_data->constant;
  }

  return NULL;
}

//
// enode is a literal if it is
// an atom or a negated atom
//
inline bool Enode::isLit( ) const
{
  if ( !isTerm( ) ) return false;
  if ( isAtom( ) ) return true;
  // if ( car->getArity( ) != ENODE_ARITY_1 ) return false;
  if ( getArity( ) != 1 ) return false;
  Enode * arg = get1st( );
  return isNot( ) && arg->isAtom( );
}
//
// enode is an atom if it has boolean type and
// it is not a boolean operator. Constants true
// and false are considered atoms
//
inline bool Enode::isAtom( ) const
{
  if ( !isTerm( ) ) return false;
  if ( !hasSortBool( ) ) return false;
  if ( isBooleanOperator( ) ) return false;

  return true;
}
//
// enode is a tatom if it has boolean type and
// it is not a boolean operator nor a boolean variable
// nor constants true and false.
//
inline bool Enode::isTAtom( ) const
{
  if ( !isTerm( ) ) return false;
  if ( !isAtom( ) ) return false;
  if ( isVar( ) )   return false;
  if ( isTrue( ) )  return false;
  if ( isFalse( ) ) return false;
  return true;
}
//
// enode is a literal if it is
// an atom or a negated atom
//
inline bool Enode::isTLit( ) const
{
  if ( !isTerm( ) ) return false;
  if ( isTAtom( ) ) return true;
  if ( getArity( ) != 1 ) return false;
  Enode * arg = get1st( );
  return isNot( ) && arg->isTAtom( );
}

inline bool Enode::isVar( ) const
{
  if ( car->getId( )    <= ENODE_ID_LAST ) return false;     // If the code is predefined, is not a variable
  if ( isConstant( ) ) return false;                         // If it's a constant is not a variable
  if ( !isTerm( ) ) return false;                            // If is not a term is not a variable
  if ( getArity( ) != 0 ) return false;
  return car->isSymb( );                                     // Final check
}

inline bool Enode::isConstant( ) const
{
  if ( !isTerm( ) ) return false;                        // Check is a term
  return isTrue( ) || isFalse( ) || car->isNumb( );      // Only numbers, true, false are constant
}

inline bool Enode::isBooleanOperator( ) const
{
  return isAnd( ) || isOr( )  || isNot( )
      || isIff( )
      || isXor( ) || isImplies( )
      ;
}

inline Enode * Enode::get1st ( ) const
{
  assert( isTerm( ) );
  assert( getArity( ) > 0 );
  return getCdr( )->getCar( );
}

inline Enode * Enode::get2nd ( ) const
{
  assert( isTerm( ) );
  assert( getArity( ) > 1 );
  return getCdr( )->getCdr( )->getCar( );
}

inline Enode * Enode::get3rd ( ) const
{
  assert( isTerm( ) );
  assert( getArity( ) > 2 );
  return getCdr( )->getCdr( )->getCdr( )->getCar( );
}

inline unsigned Enode::sizeInMem( ) const
{
  unsigned size = sizeof( Enode );
  if ( isSymb( ) )
  {
    assert( symb_data );
    assert( symb_data->name );
    size += sizeof( SymbData ) + strlen( symb_data->name );
  }
  if ( isNumb( ) )
  {
    assert( symb_data );
    assert( symb_data->name );
    assert( symb_data->value );
    size += sizeof( SymbData ) + strlen( symb_data->name ) + sizeof( double );
  }
  if ( atom_data ) size += sizeof( AtomData );
  if ( cong_data )
  {
    size += sizeof( CongData );
    if ( cong_data->term_data ) size += sizeof( TermData );
  }
  return size;
}

inline void Enode::allocCongData( )
{
  assert( isTerm( ) || isList( ) );
  assert( cong_data == NULL );
  cong_data = new CongData( id, const_cast<Enode *>(this) );
  if ( isTerm( ) )
    cong_data->term_data = new TermData( const_cast<Enode *>(this) );
}

inline void Enode::deallocCongData( )
{
  assert( isTerm( ) || isList( ) );
  assert( cong_data != NULL );
  if ( isTerm( ) )
  {
    assert( cong_data->term_data != NULL );
    delete cong_data->term_data;
    cong_data->term_data = NULL;
  }
  delete cong_data;
  cong_data = NULL;
}

inline bool Enode::isArithmeticOp( ) const
{
  //
  // Put the exhaustive list of arithmetic operators here
  //
  if ( isPlus    ( )
    || isUminus  ( )
    || isTimes   ( )
    || isLeq     ( )
    || isEq      ( ) )
    return true;

  return false;
}

inline bool Enode::isUFOp( ) const
{
  //
  // Put the exhaustive list of uf operators here
  //
  if ( isUf       ( )
    || isUp       ( )
    || isDistinct ( )
    || isEq       ( ) )
    return true;

  return false;
}
#endif
