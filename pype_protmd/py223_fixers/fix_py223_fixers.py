from lib2to3.fixer_base import BaseFix
from lib2to3.pgen2 import token


class FixPy223Fixers(BaseFix):
    '''
    '''

#    _accept_type = token.NAME

    def match(self, node):
#        print(node.value)
        if node.type==314:
#        if node.__dict__['type'] == 314:
#            print(node.__dict__)
            typs = ([x.type for x in node.children])
            if typs == [1,3]:
                return True
        return False

#         if node.value == 'raise':
# #            print(node.__dict__)
#             print(node.value)
#             try:
#                 print(node.__dict__['parent'].__dict__)
#             except:

#                 pass
#             return True
        # return False

    def transform(self, node, results):
        node.children[1].value= '''Exception("%s")'''%node.children[1].value
        node.changed()


