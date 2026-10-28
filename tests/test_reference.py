"""Scientific identity tests, including independent negative-strand and splice-junction examples."""
import sys,json,unittest
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
from reference_lib import positions,splice,blocks_for_slice,find_orfs,translate
from validate_reference import validate
class ReferenceTests(unittest.TestCase):
    def test_negative_strand_order_and_sequence(self):
        self.assertEqual(positions([[1,3],[7,9]],'-'),[9,8,7,3,2,1])
        self.assertEqual(splice('AAACCCGGG',1,[[1,3],[7,9]],'-'),'CCCTTT')
        self.assertEqual(blocks_for_slice([9,8,7,3,2,1],1,5),[[2,3],[7,8]])
    def test_nested_orfs_first_stop_and_minimum(self):
        seq='ATGATG'+'AAA'*7+'TAA'+'ATGAAATAG'
        found=list(find_orfs(seq,8));self.assertEqual(found,[(0,30),(3,30)])
        self.assertTrue(all(translate(seq[a:b]).endswith('*') for a,b in found))
    def test_current_reference(self):
        data=json.loads(Path('data/current/reference.json').read_text());loc=json.loads(Path('mycn_locus_seq.json').read_text());self.assertEqual(validate(data,loc)['status'],'PASS')
        expected={'201':'ENST00000281043.4','202':'ENST00000638417.2','203':'ENST00000703162.1','204':'ENST00000885101.2','205':'ENST00000930195.2','206':'ENST00000930196.1','207':'ENST00001138178.1','XM':'NM_001446166.1'}
        for k,v in expected.items():self.assertEqual(data['transcripts'][k]['transcript_id'],v)
        self.assertEqual(data['transcripts']['202']['exons'][0][0],15940044)
        self.assertEqual(data['transcripts']['204']['exons'][-1][1],15947022)
        self.assertEqual(data['transcripts']['205']['exons'][-1][1],15947022)
        for t in data['transcripts'].values():self.assertEqual(t['evidence_class'],'DATABASE_ANNOTATED')
class ProductTests(unittest.TestCase):
    def setUp(self):
        self.data=json.loads(Path('data/current/reference.json').read_text())
    def test_curated_identity_and_original_palette(self):
        expected={'nmyc':('201',464,'#4A6FA5'),'dmycn':('202',253,'#D98A3D'),'mycnot':('201',112,'#C0504D'),'musep':('205',86,'#B0568C')}
        import hashlib
        for key,(tx,length,color) in expected.items():
            p=self.data['products'][key];o=self.data['orfs'][p['orf_id']]
            self.assertEqual(p['primary_transcript'],tx);self.assertIn(tx,o['carriers']);self.assertEqual(p['protein_length'],length);self.assertEqual(p['color'],color);self.assertEqual(o['color'],color)
            self.assertEqual(p['protein_sequence'],translate(o['nt'])[:-1]);self.assertEqual(hashlib.sha256(p['protein_sequence'].encode()).hexdigest(),p['protein_sequence_hash'])
        self.assertIn('MNOP',self.data['products']['musep']['alternate_names'])
        ps=self.data['products'];self.assertEqual(ps['dmycn']['protein_sequence'][:52],ps['mycnot']['protein_sequence'][:52]);self.assertEqual(ps['dmycn']['protein_sequence'][:52],ps['musep']['protein_sequence'][:52]);self.assertEqual(ps['dmycn']['protein_sequence'][52:],ps['nmyc']['protein_sequence'][263:])
    def test_singleton_preservation_and_exact_window_scope(self):
        rel=self.data['relationships'];self.assertEqual(rel['counts']['orf'],{'all':42,'unique':11,'shared':31})
        for f in rel['features']:
            self.assertTrue(f['scope']);self.assertTrue(f['universe']);self.assertEqual(f['relationship']=='unique',len(f['carriers'])==1)
            if f['unit']=='protein_sequence':
                seq=self.data['products'][f['product_id']]['protein_sequence'];a,b=f['window_start_range']
                for i in range(a-1,b):
                    carriers=sorted(k for k,p in self.data['products'].items() if seq[i:i+8] in p['protein_sequence'])
                    self.assertEqual(carriers,f['carriers'])
class UpdateTests(unittest.TestCase):
    def test_release_and_refseq_detection(self):
        from check_for_reference_updates import check
        page='<h1>Release 51 (GRCh38.p14)</h1>'
        from Bio.SeqRecord import SeqRecord
        from Bio.Seq import Seq
        from Bio import SeqIO
        from io import StringIO
        fixture=SeqRecord(Seq('ATG'),id='NM_001446166.1',name='NM_001446166',description='Test fixture.')
        fixture.annotations={'molecule_type':'mRNA','comment':'REVIEWED REFSEQ'}
        buffer=StringIO();SeqIO.write(fixture,buffer,'genbank');record=buffer.getvalue()
        manifest={'gencode_release':50,'refseq_accession':'NM_001446166.1'}
        result=check(manifest,page,record);self.assertEqual(result['status'],'UPDATE_DETECTED');self.assertFalse(result['automatic_publication'])
        self.assertEqual(check(manifest,page.replace('51','50'),record)['status'],'NO_UPDATE')
if __name__=='__main__':unittest.main()
