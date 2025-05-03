Attribute VB_Name = "Module2"
Sub Converter_selection()
Attribute Converter_selection.VB_ProcData.VB_Invoke_Func = " \n14"
    'MsgBox (" O programa esta a iniciar ")
    'Application.MacroOptions Macro:="Converter_selection", Description:="", _
        'ShortcutKey:=""
    'Application.Goto Reference:="Converter_selection"
    'ActiveSheet.Shapes.Range(Array("CommandButton1")).Select
    'Selection.Verb Verb:=xlPrimary
    'Range("B3").Select
    'ActiveSheet.Shapes.Range(Array("CommandButton1")).Select
    UserForm1.Show
    ActiveWorkbook.Save
    Application.WindowState = xlNormal
End Sub
