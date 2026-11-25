"""Portable patent review table with direct links to standalone local JSON."""
import argparse
import csv
import html
from pathlib import Path
from urllib.parse import quote

def main():
    p=argparse.ArgumentParser();p.add_argument('--results',type=Path,required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args()
    with (a.results/'patents.csv').open(encoding='utf-8-sig') as f:
        rows=[r for r in csv.DictReader(f) if r['status']=='processed']
    headers=['专利号','标题 / 公司','书目与引证JSON','申请信息 / 发明人','摘要与其他元数据']
    visible={'Patent','title','company_name','bibliography_json_path','application_number','application_date','inventor_names','abstract','status'}
    labels={'abstract_language':'摘要语言','is_pct':'是否识别到PCT关联','pct_number':'PCT申请号','pct_relationship':'PCT与本申请关系','wo_publication_number':'WO国际公开号','wo_publication_date':'WO国际公开日','priority_date':'最早记载优先权日','priority_date_status':'优先权日期状态','application_publication_number':'美国本申请公开号','application_publication_date':'美国本申请公开日','grant_date':'授权/再授权日','grant_date_type':'授权日期类型','indication':'适应症/专利披露治疗用途','other_disclosed_indications':'其他层级披露用途','ipc_codes':'IPC','cpc_codes':'CPC','patent_category':'专利内容分类','technical_field':'技术领域','applicant_all_names':'申请人','inventor_source':'发明人角色来源','parent_application_numbers':'关联母案申请号','parent_application_relationships':'母案关联类型','patent_citation_count':'专利引证条数','npl_citation_count':'非专利引证条数','document_metadata_review_reason':'全文字段待核对项','archive':'来源ZIP','xml_file':'来源XML'}
    def esc(x):return html.escape(str(x or ''))
    def link(r,key,text):
        relative=r[key];destination=a.results/relative
        assert destination.is_file() and destination.resolve().is_relative_to(a.results.resolve())
        return '<a target="_blank" rel="noopener" href="'+quote(relative,safe='/')+'">'+text+'</a><small>'+esc(relative)+'</small>'
    body=[]
    for r in rows:
        extras=''.join('<dt>'+esc(labels.get(k,k))+'</dt><dd>'+esc(v)+'</dd>' for k,v in r.items() if k not in visible and v)
        body.append('<tr><td class="patent">'+esc(r['Patent'])+'</td><td><strong>'+esc(r['title'])+'</strong><small>'+esc(r['company_name'])+'</small></td><td>'+link(r,'bibliography_json_path','打开书目及引证 JSON')+'</td><td>'+esc(r['application_number'])+'<br>'+esc(r['application_date'])+'<small>'+esc(r['inventor_names'])+'</small></td><td><details><summary>查看摘要</summary><p>'+esc(r['abstract'])+'</p></details><details><summary>查看其他元数据</summary><dl>'+extras+'</dl></details></td></tr>')
    page='''<!doctype html><html lang="zh-CN"><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>USPTO专利样例审阅表</title><style>
body{margin:24px;background:#f6f8fa;color:#202833;font:14px/1.55 system-ui,-apple-system,sans-serif}h1{font-size:25px;color:#244663;margin-bottom:8px}.intro{max-width:1100px}.wrap{max-height:78vh;overflow:auto;background:white;border:1px solid #d6e0e9;border-radius:8px}table{border-collapse:separate;border-spacing:0;width:100%;min-width:1350px}th{position:sticky;top:0;background:#244663;color:white;text-align:left;z-index:2;padding:14px 12px}td{padding:14px 12px;vertical-align:top;border-bottom:1px solid #dae2e9}tr:nth-child(even){background:#f4f8fb}td:nth-child(2){min-width:290px;max-width:400px}td:nth-child(3){min-width:230px}td:last-child{min-width:260px;max-width:480px}.patent{white-space:nowrap;font-weight:600}small{display:block;color:#5d6a79;margin-top:8px;overflow-wrap:anywhere}a{color:#1263a3}summary{cursor:pointer;color:#1263a3;padding:4px 0}p,dd{overflow-wrap:anywhere}dt{font-weight:600;margin-top:10px}dd{margin-left:0;font-size:13px}footer{color:#596674;margin-top:18px;max-width:1200px}
</style><h1>USPTO 医药专利样例</h1><p class="intro">'''+str(len(rows))+''' 件医药专利 · 在 open 服务器以 20 进程运行。摘要、当事人、关联申请和引证明细逐专利保存为 JSON。请完整解压，保持本表与 bibliography 目录同级；在桌面浏览器中点击蓝色链接查看文件。</p><div class="wrap"><table><thead><tr>'''+''.join('<th>'+h+'</th>' for h in headers)+'''</tr></thead><tbody>'''+''.join(body)+'''</tbody></table></div><footer>当事人及关联申请来自本件文献记载；引证不作为本申请日期或母案信息。当前为40件医药样例，全量未运行。</footer></html>'''
    a.output.parent.mkdir(parents=True,exist_ok=True);a.output.write_text(page,encoding='utf-8')
    print(f'HTML review: {len(rows)} patents, {len(rows)} local JSON links, {a.output}')

if __name__=='__main__':main()
