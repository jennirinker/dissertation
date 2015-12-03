"""
testing fast 7 dictionary creation
"""
import jr_fast

fast_fpath = 'C:\\Users\\jrinker\\Documents\\GitHub\\' + \
                'dissertation\\FAST_models\\FAST7\\' + \
                'WP0.75A08V00\\WP0.75A08V00_24134.fst'

turb_dict = jr_fast.CreateFAST7Dict(fast_fpath)


