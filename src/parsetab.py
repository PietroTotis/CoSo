
# parsetab.py
# This file is automatically generated. Do not edit.
# pylint: disable=W,C,R
_tabversion = '3.10'

_lr_method = 'LALR'

_lr_signature = 'COMMA COUNT DIFF EQUALS GT IN INTER LAB LABEL LPAR LRPAR LSPAR LT MAX MIN NOT NUMBER PART PROP REPEAT RPAR RRPAR RSPAR SEMI SUM UNION UNIVERSEprogram : statement\n        | statement program\n        statement : declare_set SEMI\n        | arrangement SEMI\n        | aggcmp SEMI\n        | pos_constraint SEMI\n        | count_constraint SEMI\n        entity : NUMBER\n        | LABEL\n        entity_list : entity\n        | entity COMMA entity_list\n        comp : EQUALS\n        | LT\n        | GT\n        | GT EQUALS\n        | LT EQUALS\n        | DIFF EQUALS\n        base_set : LABEL\n        | LABEL LSPAR NUMBER RSPAR\n        | UNIVERSE\n        set : base_set\n        | PART\n        | LRPAR set RRPAR\n        | NOT set\n        | set INTER set\n        | set UNION set\n        declare_set : PROP LABEL EQUALS LPAR entity_list RPAR\n        | UNIVERSE LABEL EQUALS LPAR entity_list RPAR\n        | PROP UNIVERSE EQUALS LPAR entity_list RPAR\n        | LAB PROP LABEL EQUALS LPAR entity_list RPAR\n        | PROP LABEL\n        | LAB PROP LABEL\n        arrangement : LABEL IN LPAR REPEAT set RPAR\n        | LABEL IN LPAR set RPAR\n        | LABEL IN LSPAR REPEAT set RSPAR\n        | LABEL IN LSPAR set RSPAR\n        | LABEL IN LPAR LPAR set RPAR RPAR\n        | LABEL IN LSPAR LPAR set RPAR RSPAR\n        math_op : SUM\n        | MIN\n        | MAX\n        aggcmp  : math_op LRPAR set RRPAR comp NUMBERpos_constraint : LABEL LSPAR NUMBER RSPAR EQUALS entity\n        | LABEL LSPAR NUMBER RSPAR IN set\n        | LABEL LSPAR NUMBER RSPAR EQUALS set\n        | LABEL LSPAR NUMBER RSPAR EQUALS LPAR entity_list RPAR\n        count_constraint : COUNT set comp NUMBER\n        | COUNT LRPAR count_constraint RRPAR comp NUMBER\n        '
    
_lr_action_items = {'PROP':([0,2,11,18,19,20,21,22,],[8,8,28,-3,-4,-5,-6,-7,]),'UNIVERSE':([0,2,8,13,18,19,20,21,22,29,31,34,39,40,44,47,48,59,60,62,64,89,90,],[10,10,24,36,-3,-4,-5,-6,-7,36,36,36,36,36,36,36,36,36,36,36,36,36,36,]),'LAB':([0,2,18,19,20,21,22,],[11,11,-3,-4,-5,-6,-7,]),'LABEL':([0,2,8,10,13,18,19,20,21,22,28,29,31,34,39,40,44,47,48,57,58,59,60,62,64,66,89,90,92,97,106,],[9,9,23,27,35,-3,-4,-5,-6,-7,43,35,35,35,35,35,35,35,35,78,78,35,35,35,35,78,103,35,78,78,78,]),'COUNT':([0,2,18,19,20,21,22,31,],[13,13,-3,-4,-5,-6,-7,13,]),'SUM':([0,2,18,19,20,21,22,],[14,14,-3,-4,-5,-6,-7,]),'MIN':([0,2,18,19,20,21,22,],[15,15,-3,-4,-5,-6,-7,]),'MAX':([0,2,18,19,20,21,22,],[16,16,-3,-4,-5,-6,-7,]),'$end':([1,2,17,18,19,20,21,22,],[0,-1,-2,-3,-4,-5,-6,-7,]),'SEMI':([3,4,5,6,7,23,32,33,35,36,43,55,69,70,71,76,81,85,87,95,96,98,100,101,103,104,105,107,108,110,111,113,114,116,117,],[18,19,20,21,22,-31,-21,-22,-18,-20,-32,-24,-47,-25,-26,-23,-8,-34,-36,-19,-27,-29,-33,-35,-9,-43,-45,-44,-28,-42,-48,-37,-38,-30,-46,]),'IN':([9,65,],[25,90,]),'LSPAR':([9,25,35,103,],[26,40,56,56,]),'LRPAR':([12,13,14,15,16,29,31,34,39,40,44,47,48,59,60,62,64,89,90,],[29,31,-39,-40,-41,44,44,44,44,44,44,44,44,44,44,44,44,44,44,]),'PART':([13,29,31,34,39,40,44,47,48,59,60,62,64,89,90,],[33,33,33,33,33,33,33,33,33,33,33,33,33,33,33,]),'NOT':([13,29,31,34,39,40,44,47,48,59,60,62,64,89,90,],[34,34,34,34,34,34,34,34,34,34,34,34,34,34,34,]),'EQUALS':([23,24,27,30,32,33,35,36,43,50,51,52,55,65,68,70,71,75,76,95,],[37,38,42,49,-21,-22,-18,-20,67,72,73,74,-24,89,49,-25,-26,49,-23,-19,]),'LPAR':([25,37,38,39,40,42,67,89,],[39,57,58,59,64,66,92,106,]),'NUMBER':([26,46,49,50,51,56,57,58,66,72,73,74,89,92,93,94,97,106,],[41,69,-12,-13,-14,77,81,81,81,-16,-15,-17,81,81,110,111,81,81,]),'INTER':([30,32,33,35,36,45,54,55,61,63,70,71,76,83,84,86,88,95,103,105,107,],[47,-21,-22,-18,-20,47,47,47,47,47,47,47,-23,47,47,47,47,-19,-18,47,47,]),'UNION':([30,32,33,35,36,45,54,55,61,63,70,71,76,83,84,86,88,95,103,105,107,],[48,-21,-22,-18,-20,48,48,48,48,48,48,48,-23,48,48,48,48,-19,-18,48,48,]),'LT':([30,32,33,35,36,55,68,70,71,75,76,95,],[50,-21,-22,-18,-20,-24,50,-25,-26,50,-23,-19,]),'GT':([30,32,33,35,36,55,68,70,71,75,76,95,],[51,-21,-22,-18,-20,-24,51,-25,-26,51,-23,-19,]),'DIFF':([30,32,33,35,36,55,68,70,71,75,76,95,],[52,-21,-22,-18,-20,-24,52,-25,-26,52,-23,-19,]),'RRPAR':([32,33,35,36,45,53,54,55,69,70,71,76,95,111,],[-21,-22,-18,-20,68,75,76,-24,-47,-25,-26,-23,-19,-48,]),'RPAR':([32,33,35,36,55,61,70,71,76,78,79,80,81,82,83,84,88,91,95,99,109,112,115,],[-21,-22,-18,-20,-24,85,-25,-26,-23,-9,96,-10,-8,98,99,100,102,108,-19,113,116,-11,117,]),'RSPAR':([32,33,35,36,41,55,63,70,71,76,77,86,95,102,],[-21,-22,-18,-20,65,-24,87,-25,-26,-23,95,101,-19,114,]),'REPEAT':([39,40,],[60,62,]),'COMMA':([78,80,81,],[-9,97,-8,]),}

_lr_action = {}
for _k, _v in _lr_action_items.items():
   for _x,_y in zip(_v[0],_v[1]):
      if not _x in _lr_action:  _lr_action[_x] = {}
      _lr_action[_x][_k] = _y
del _lr_action_items

_lr_goto_items = {'program':([0,2,],[1,17,]),'statement':([0,2,],[2,2,]),'declare_set':([0,2,],[3,3,]),'arrangement':([0,2,],[4,4,]),'aggcmp':([0,2,],[5,5,]),'pos_constraint':([0,2,],[6,6,]),'count_constraint':([0,2,31,],[7,7,53,]),'math_op':([0,2,],[12,12,]),'set':([13,29,31,34,39,40,44,47,48,59,60,62,64,89,90,],[30,45,54,55,61,63,54,70,71,83,84,86,88,105,107,]),'base_set':([13,29,31,34,39,40,44,47,48,59,60,62,64,89,90,],[32,32,32,32,32,32,32,32,32,32,32,32,32,32,32,]),'comp':([30,68,75,],[46,93,94,]),'entity_list':([57,58,66,92,97,106,],[79,82,91,109,112,115,]),'entity':([57,58,66,89,92,97,106,],[80,80,80,104,80,80,80,]),}

_lr_goto = {}
for _k, _v in _lr_goto_items.items():
   for _x, _y in zip(_v[0], _v[1]):
       if not _x in _lr_goto: _lr_goto[_x] = {}
       _lr_goto[_x][_k] = _y
del _lr_goto_items
_lr_productions = [
  ("S' -> program","S'",1,None,None,None),
  ('program -> statement','program',1,'p_program','parser.py',121),
  ('program -> statement program','program',2,'p_program','parser.py',122),
  ('statement -> declare_set SEMI','statement',2,'p_statement','parser.py',126),
  ('statement -> arrangement SEMI','statement',2,'p_statement','parser.py',127),
  ('statement -> aggcmp SEMI','statement',2,'p_statement','parser.py',128),
  ('statement -> pos_constraint SEMI','statement',2,'p_statement','parser.py',129),
  ('statement -> count_constraint SEMI','statement',2,'p_statement','parser.py',130),
  ('entity -> NUMBER','entity',1,'p_entity','parser.py',134),
  ('entity -> LABEL','entity',1,'p_entity','parser.py',135),
  ('entity_list -> entity','entity_list',1,'p_entity_list','parser.py',146),
  ('entity_list -> entity COMMA entity_list','entity_list',3,'p_entity_list','parser.py',147),
  ('comp -> EQUALS','comp',1,'p_comp','parser.py',155),
  ('comp -> LT','comp',1,'p_comp','parser.py',156),
  ('comp -> GT','comp',1,'p_comp','parser.py',157),
  ('comp -> GT EQUALS','comp',2,'p_comp','parser.py',158),
  ('comp -> LT EQUALS','comp',2,'p_comp','parser.py',159),
  ('comp -> DIFF EQUALS','comp',2,'p_comp','parser.py',160),
  ('base_set -> LABEL','base_set',1,'p_base_set','parser.py',168),
  ('base_set -> LABEL LSPAR NUMBER RSPAR','base_set',4,'p_base_set','parser.py',169),
  ('base_set -> UNIVERSE','base_set',1,'p_base_set','parser.py',170),
  ('set -> base_set','set',1,'p_set','parser.py',179),
  ('set -> PART','set',1,'p_set','parser.py',180),
  ('set -> LRPAR set RRPAR','set',3,'p_set','parser.py',181),
  ('set -> NOT set','set',2,'p_set','parser.py',182),
  ('set -> set INTER set','set',3,'p_set','parser.py',183),
  ('set -> set UNION set','set',3,'p_set','parser.py',184),
  ('declare_set -> PROP LABEL EQUALS LPAR entity_list RPAR','declare_set',6,'p_declare_set','parser.py',214),
  ('declare_set -> UNIVERSE LABEL EQUALS LPAR entity_list RPAR','declare_set',6,'p_declare_set','parser.py',215),
  ('declare_set -> PROP UNIVERSE EQUALS LPAR entity_list RPAR','declare_set',6,'p_declare_set','parser.py',216),
  ('declare_set -> LAB PROP LABEL EQUALS LPAR entity_list RPAR','declare_set',7,'p_declare_set','parser.py',217),
  ('declare_set -> PROP LABEL','declare_set',2,'p_declare_set','parser.py',218),
  ('declare_set -> LAB PROP LABEL','declare_set',3,'p_declare_set','parser.py',219),
  ('arrangement -> LABEL IN LPAR REPEAT set RPAR','arrangement',6,'p_arrangement','parser.py',243),
  ('arrangement -> LABEL IN LPAR set RPAR','arrangement',5,'p_arrangement','parser.py',244),
  ('arrangement -> LABEL IN LSPAR REPEAT set RSPAR','arrangement',6,'p_arrangement','parser.py',245),
  ('arrangement -> LABEL IN LSPAR set RSPAR','arrangement',5,'p_arrangement','parser.py',246),
  ('arrangement -> LABEL IN LPAR LPAR set RPAR RPAR','arrangement',7,'p_arrangement','parser.py',247),
  ('arrangement -> LABEL IN LSPAR LPAR set RPAR RSPAR','arrangement',7,'p_arrangement','parser.py',248),
  ('math_op -> SUM','math_op',1,'p_math_op','parser.py',291),
  ('math_op -> MIN','math_op',1,'p_math_op','parser.py',292),
  ('math_op -> MAX','math_op',1,'p_math_op','parser.py',293),
  ('aggcmp -> math_op LRPAR set RRPAR comp NUMBER','aggcmp',6,'p_aggcmp','parser.py',298),
  ('pos_constraint -> LABEL LSPAR NUMBER RSPAR EQUALS entity','pos_constraint',6,'p_pos_constraint','parser.py',312),
  ('pos_constraint -> LABEL LSPAR NUMBER RSPAR IN set','pos_constraint',6,'p_pos_constraint','parser.py',313),
  ('pos_constraint -> LABEL LSPAR NUMBER RSPAR EQUALS set','pos_constraint',6,'p_pos_constraint','parser.py',314),
  ('pos_constraint -> LABEL LSPAR NUMBER RSPAR EQUALS LPAR entity_list RPAR','pos_constraint',8,'p_pos_constraint','parser.py',315),
  ('count_constraint -> COUNT set comp NUMBER','count_constraint',4,'p_count_constraint','parser.py',344),
  ('count_constraint -> COUNT LRPAR count_constraint RRPAR comp NUMBER','count_constraint',6,'p_count_constraint','parser.py',345),
]
