import pandas as pd
import click


def check_callers(df):
    list_callers = list()
    if df.loc[0, 'SAGE'] == True:
        list_callers.extend(['SAGE'])
    if df.loc[0, 'STRELKA'] == True:
        list_callers.extend(['STRELKA'])
    if df.loc[0, 'MUTECT'] == True:
        list_callers.extend(['MUTECT'])
    df['Callers_intersection'] = ",".join(list_callers)
    return (df)


@click.command()

@click.option('--maf_sage',
              '-sa',
              required = True,
              help="Input must be the path to the processed MAF file from SAGE calling")
@click.option('--maf_mutect',
              '-mu',
              required = True,
              help="Input must be the path to the processed MAF file from MuTect2 calling")
@click.option('--maf_strelka',
              '-st',
              required = True,
              help="Input must be the path to the processed MAF file from Strelka calling")
@click.option('--output_dire',
              '-o',
              required = True,
              help="Output MAF with variants in 2 out of the 3 callers")
@click.option('--sample_name',
              '-sn',
              required = True,
              help="Sample name to put in SAMPLE column. This will be as well the name of the file")
@click.option('--by_chrom/--single_file',
              '-c/-s',
              required = False,
              default = False,
              show_default=True,
              help="Whether the output should be in different files by chromosome. Default is False")


def cli(maf_sage, maf_mutect, maf_strelka, output_dire, sample_name, by_chrom):


    # READ RESULTS
    df_maf_sage = pd.read_csv(maf_sage, sep='\t',compression='gzip')
    df_maf_mutect = pd.read_csv(maf_mutect, sep='\t',compression='gzip')
    df_maf_strelka = pd.read_csv(maf_strelka, sep='\t',compression='gzip')

    # CREATE COLUMN OF VARIANT FOR INTERSECTION
    df_maf_sage['Variant'] = df_maf_sage.apply(lambda x: x['#CHROM']+str(x['POS'])+x['REF']+x['ALT'], axis=1)
    df_maf_mutect['Variant'] = df_maf_mutect.apply(lambda x: x['#CHROM'] + str(x['POS']) + x['REF'] + x['ALT'], axis=1)
    df_maf_strelka['Variant'] = df_maf_strelka.apply(lambda x: x['#CHROM'] + str(x['POS']) + x['REF'] + x['ALT'], axis=1)

    # INTERSECT SETS
    sage_mutect = set(df_maf_sage['Variant'].tolist()).intersection(set(df_maf_mutect['Variant'].tolist()))
    sage_strelka = set(df_maf_sage['Variant'].tolist()).intersection(set(df_maf_strelka['Variant'].tolist()))
    strelka_mutect = set(df_maf_strelka['Variant'].tolist()).intersection(set(df_maf_mutect['Variant'].tolist()))

    # RECOVER MUTATIONS AND DECIDE COLUMNS TO KEEP (SAGE>STRELKA>MUTECT2)
    df1 = df_maf_sage[df_maf_sage['Variant'].isin(sage_mutect)]
    df1['SAGE'] = True
    df1['MUTECT'] = True

    df2 = df_maf_sage[df_maf_sage['Variant'].isin(sage_strelka)]
    df2['SAGE'] = True
    df2['STRELKA'] = True

    df3 = df_maf_strelka[df_maf_strelka['Variant'].isin(strelka_mutect)]
    df3['STRELKA'] = True
    df3['MUTECT'] = True

    df = df1.append(df2, ignore_index=True, sort=False)
    df = df.append(df3, ignore_index=True, sort=False)

    df['SAGE'].fillna(False, inplace=True)
    df['STRELKA'].fillna(False, inplace=True)
    df['MUTECT'].fillna(False, inplace=True)

    # REMOVE REPEATED
    grps = df.groupby('Variant')

    dff = pd.DataFrame()

    for g in grps.groups:
        df_var = grps.get_group(g).reset_index()
        if len(df_var) == 1:
            df_var = check_callers(df_var)
            dff = dff.append(df_var, ignore_index=True, sort=False)
        elif len(df_var) > 1:
            if (len(df_var) == 3) and (True in df_var['SAGE'].tolist()):
                if True in df_var['SAGE'].tolist():
                    df_var = df_var[df_var['SAGE'] == True]
                    df_var.drop_duplicates( subset=['SAGE'], keep='first', inplace=True)
                    df_var['Callers_intersection'] = 'SAGE,STRELKA,MUTECT'
                    dff = dff.append(df_var, ignore_index=True, sort=False)
                else:
                    print("equal to 3 without SAGE ?. Impossible")
                    print(df_var)
            elif len(df_var) == 2:
                if (True in df_var['SAGE'].tolist()) and (True in df_var['STRELKA'].tolist()):
                    df_var = df_var[df_var['SAGE'] == True]
                    df_var['Callers_intersection'] = 'SAGE,STRELKA'
                elif (True in df_var['SAGE'].tolist()) and (True in df_var['MUTECT'].tolist()):
                    df_var = df_var[df_var['SAGE'] == True]
                    df_var['Callers_intersection'] = 'SAGE,MUTECT'
                elif (not True in df_var['SAGE'].tolist()) and (True in df_var['STRELKA'].tolist()):
                    df_var = df_var[df_var['STRELKA'] == True]
                    df_var['Callers_intersection'] = 'STRELKA,MUTECT'
                else:
                    print("another scenario ? show me")
                    print(df_var)
                dff = dff.append(df_var, ignore_index=True, sort=False)
            else:
                print("This should be printed")
        else:
            pass
    
    # DROP INDEX, SAGE, STRELKA and MUTECT COLUMNS
    dff.drop(columns=['index','SAGE', 'STRELKA', 'MUTECT'], inplace=True)
    
    # ADD SAMPLE COLUMN
    dff['SAMPLE'] = sample_name
    
    if by_chrom == True:
    
        # SAVE TABLES DIVIDED BY CHROMOSOME
        chrom_groups = dff.groupby("#CHROM")

        for chr_ in chrom_groups.groups:
            df_chr = chrom_groups.get_group(chr_)            
            df_chr.to_csv(output_dire+sample_name+"_"+chr_+'.maf.gz', sep='\t', index=False,compression='gzip')
            
    elif by_chrom == False:

        # SAVE ONE TABLE
        dff.to_csv(output_dire+sample_name+'.maf.gz', sep='\t', index=False,compression='gzip')
        
    else:
        print('by_chrom should be a True or False statement')


if __name__ == '__main__':
    cli()